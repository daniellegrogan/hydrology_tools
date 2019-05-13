#!/usr/bin/perl -w

#######################################################################
#
#	This code generates a Cell Table for the Water Balance Model (WBM).
#
#	Written by Dr. A. Prusevich (alex.proussevitch@unh.edu)
#
#	January 2011
#	Modified-	February 2019
#
#######################################################################

use strict;
use Benchmark;
use Clone qw/clone/;
use Getopt::Long;
use Geo::GDAL;
use File::Basename;
use File::Temp;
use FileHandle;
use Math::Trig qw/pi/;
use Math::VecStat;
use PDL;
use PDL::IO::FlexRaw;
use PDL::NiceSlice;
use PDL::NetCDF;
use Fcntl;
use Inline qw/Pdlpp/;
use RIMS;

my ($init_file, $io_file)	= (script_dir().'wbm_path.init', script_dir().'wbm_io.pl');
{			### wbm_io.pl must be in the same directory
  local @ARGV	= ($init_file);
  require $io_file;
}			### wbm_path.init must be in the same directory

			### Use BIGPDL if the network is really big (over 1GB)
$PDL::BIGPDL = 1 unless $PDL::BIGPDL;

use vars qw(*OLDERR);		# To avoid silly message-
open OLDERR, ">&STDERR";	# "No UNIDATA NC_GLOBAL:Conventions attribute"
STDOUT->autoflush(1);				# Disable buffering

my $time_start = time();
#######################################################################
#############   Process and check command line inputs   ###############

my ($help,$verbose,$subset,$runInt,$file_ID,$proj) = (0,0,-1,0,'','epsg:4326');
						# Get command line options
usage() if !GetOptions('h'=>\$help, 'v'=>\$verbose, 'ib'=>\$runInt, 'b=s'=>\$file_ID,
	'p=s'=>\$proj, 'sub=i'=>\$subset) or $help;

 my $file_ntwk = shift() or usage();
 my $file_ctbl = shift();
(   $file_ctbl = $file_ntwk) =~ s/\.\w+$/.dat/ unless $file_ctbl;
(my $file_ctab = $file_ctbl) =~ s/\.\w+$/.csv/;
(my $file_upAr = $file_ctbl) =~ s/\.\w+$/_upstrArea.asc/;
(my $file_IDs  = $file_ctbl) =~ s/\.\w+$/_IDs.asc/;
(my $file_Int  = $file_ctbl) =~ s/\.\w+$/_Int.asc/;	# IDs  of internal basins

	### File names for the Subset Network
(my $file_Net		= $file_ctbl) =~ s/\.\w+$/_Subset.asc/;			# FlowDir
(my $file_ctNet		= $file_ctbl) =~ s/\.\w+$/_Subset.csv/;			# Cell table
(my $file_ctbNet	= $file_ctbl) =~ s/\.\w+$/_Subset.dat/;
(my $file_bID		= $file_ctbl) =~ s/\.\w+$/_Subset_IDs.asc/;		# Basin IDs
(my $file_upArNet	= $file_ctbl) =~ s/\.\w+$/_Subset_upstrArea.asc/;	# Upstream Area

#######################################################################
#############   Basin subsetting parameters   #########################
			### Choice for Subset # 0
my $ID		= pdl(20, 26);		# List of basin IDs
			### Choice for Subset # 1 (single point if Min=Max)
my %region      = ('lonMin'     =>  -80.05,     'lonMax'        => -78.45,
                   'latMin'     =>  44.838,     'latMax'        =>  45.782);
# my %region	= ('lonMin'	=> -71.44,	'lonMax'	=>  -70.85,
# 		   'latMin'	=>   42.94,	'latMax'	=>  43.31);
# my %region	= ('lonMin'	=> -114.5,	'lonMax'	=> -109.0,	# Arizona
# 		   'latMin'	=>   31.0,	'latMax'	=>   37.0);
# my %region	= ('lonMin'	=> -105,	'lonMax'	=> -102,	# Churchill + Nelson
# 		   'latMin'	=> 52,	'latMax'	=> 56);
# my %region 	= ('lonMin'	=> -20,		'lonMax'	=> 25,		# Wascal climate box
# 		   'latMin'	=>  0,		'latMax'	=>  30);
# my %region	= ('lonMin'	=> 20,	'lonMax'	=> 21,
# 		   'latMin'	=> 45,	'latMax'	=> 46);
			### Choice for Subset # 2
my $polygonFile	= '/net/nfs/yukon/raid5/userdata/dwisser/Networks/ecowas_5min.tif';
my $polygonID	= pdl(384);	# List of polygon IDs

#    $subset	= 1;	#<0- no subset;   0- list of IDs; 1- region; 2- polygon
my $rule	= 1;	# 0- intersect;   1- inside region/polygon
my $trim	= 2;	# 0- do not trim; 1- trim by region;         2- trim nodata
   $trim	= 0 if $subset < 0;

#######################################################################
#############   Actual Work is here   #################################

print "Building Cell Table-\n" if $verbose;

my $extent		= get_extent( $file_ntwk,$proj);
my $flowDir		= $$extent{mask};
my($lon,$lat)		= lonLat($extent);
my($colTo,$rowTo)	= flow_to($flowDir);		# Modifies $flowDir too
my $basinID		= (-e $file_ID) ? read_raster($file_ID,1,[1,0]) : 0;
my @time		= (Benchmark->new);

		### Check curcularity of the Network
printf "\tChecking curcularity of the Network (%d cells)...\n",$flowDir->ngood if $verbose;
check_circularityC($flowDir, $$extent{gTransform});			# C-coded version
# check_circularity($flowDir, $colTo, $rowTo, $$extent{gTransform});	# Plain Perl version
push @time, Benchmark->new;
printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;

		### Calculate upstream area
print "\tCalculating upstream area...\n" if $verbose;
my $up_area	= byte($flowDir)->upstrAccumAll(cell_area($lon,$lat,$extent));
push @time, Benchmark->new;
printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;

		### Build Cell table
print "\tShaping Cell table...\n" if $verbose;
my ($area,$cell_pdl) = cell_table($up_area,$flowDir,$colTo,$rowTo);
push @time, Benchmark->new;
printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;

		### Build basin IDs
print "\tBuilding basin IDs...\n" if $verbose;
   $basinID	= byte($flowDir)->basinIDs($cell_pdl) unless ref($basinID);
push @time, Benchmark->new;
printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;

	###########################################
	###   Subsetting basin network
	###########################################
my $mask= zeroes(byte,$flowDir->dims); my $IDout;
my @box	= (colRow($extent,$region{lonMin},$region{latMin}),
	   colRow($extent,$region{lonMax},$region{latMax}));
my @dim	= (0, $mask->dim(0)-1, 0, $mask->dim(1)-1);

if ($subset == 1) {	# Mask network by region
			# But check subsetting extents first
  die "\nThe subsetting regions is outside of the Network domain. Aborting...\n"
  if  $box[0] > $dim[1] || $box[2] < $dim[0] || $box[3] > $dim[3] || $box[1] < $dim[2];
  if ($box[0] < $dim[0]) { $box[0] = $dim[0]; print "\tSubsetting region min Longitute is clipped.\n" if $verbose; }
  if ($box[2] > $dim[1]) { $box[2] = $dim[1]; print "\tSubsetting region min Longitute is clipped.\n" if $verbose; }
  if ($box[3] < $dim[2]) { $box[3] = $dim[2]; print "\tSubsetting region max Latitude  is clipped.\n" if $verbose; }
  if ($box[1] > $dim[3]) { $box[2] = $dim[1]; print "\tSubsetting region max Latitude  is clipped.\n" if $verbose; }

  $mask($box[0]:$box[2],$box[3]:$box[1]) .= 1;
}
if ($subset == 2) {	# Mask network by polygon
  $mask = read_GDAL($extent,
	{'Var_Scale'=>1,'Var_Offset'=>0,'Processing'=>'','Projection'=>'epsg:4326'},
	0,$polygonFile,1,-9999)->in($polygonID);
}
if ($subset < 0) {
  $ID	= undef;
}			# Find basin IDs in mask region/polygon
elsif ($subset > 0) {
  $ID	= condition($mask, $basinID, -9999)->setvaltobad(-9999)->uniq;
  $IDout= condition($mask, -9999, $basinID)->setvaltobad(-9999)->uniq;

  $ID	= setops($ID,'XOR',setops($ID,'AND',$IDout)) if $rule;
}
			# Subset Network by basin IDs
delete $$extent{mask};		# clone does not work with PDL objects
my $extNet;
my $basinNet	= condition($basinID->in($ID), $flowDir, -9999)->setvaltobad(-9999);
my $bsnIDNet	= condition($basinID->in($ID), $basinID, -9999)->setvaltobad(-9999);

($bsnIDNet,$extNet) = trim($bsnIDNet, clone($extent), \@box, $trim);
($basinNet,$extNet) = trim($basinNet, clone($extent), \@box, $trim);

	###########################################
	###   Build Cell table for the Subset Network
	###########################################
if ($subset >= 0) {
		### Calculate upstream area
  my $extent		= get_extent( $file_ntwk,$proj);
  my $flowDir		= $$extent{mask};
  my($lon,$lat)		= lonLat($extent);
  my($colTo,$rowTo)	= flow_to($flowDir);		# Modifies $flowDir too

  print "\tCalculating upstream area for the Subset Network...\n" if $verbose;
  my $up_areaNet	= byte($flowDir)->upstrAccumAll(cell_area($lon,$lat,$extent));
  push @time, Benchmark->new;
  printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;

		### Build Cell table
  print "\tShaping Cell table for the Subset Network...\n" if $verbose;
  my ($area,$cell_pdl) = cell_table($up_area,$flowDir,$colTo,$rowTo);
  push @time, Benchmark->new;
  printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;
}
	###########################################
	###   Locate internal basins
	###########################################
if ($runInt) {
  print "\tBuilding internal basin IDs..." if $verbose;
  my $basinInt    = internal_basins($basinID, $flowDir);
  write_gridascii($file_Int, $basinInt, $extent, {FORMAT => '%d'});
  printf " %.2f %% of the Network are internal basins.\n",
	100*($basinInt->ngood/$basinID->ngood) if $verbose;
  push @time, Benchmark->new;
  printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;
}

#######################################################################

		### Build Cell table for CSV
if ($subset < 0) {
  my @area	= $area->list;
  my $cell_table= repack_data([$cell_pdl->list], [scalar(@area),5]);
  map push(@{$$cell_table[$_]}, Geo::GDAL::ApplyGeoTransform($$extent{gTransform},
	     $$cell_table[$_][3]+0.5, $$cell_table[$_][4]+0.5), $area[$_]), 0..$#$cell_table;

		### Save data to files
  write_cellTable($file_ctab,['ToCell','ToCellX','ToCellY','X','Y','Lon','Lat','SubbasinArea'],$cell_table);
  write_gridascii($file_upAr,$up_area, $extent, {FORMAT => '%.1f'});
  write_gridascii($file_IDs, $basinID, $extent, {FORMAT => '%d'});
  writeflex(	  $file_ctbl,$up_area, $area, $cell_pdl);
}
		### Save subset data to files
else {
  write_gridascii($file_Net, $basinNet,$extNet, {FORMAT => '%d'});

		### Build Cell table
  print "\tBuilding Cell table for the Subset Network...\n" if $verbose;
  my $extent		= get_extent( $file_Net,$proj);
  my $flowDir		= $$extent{mask};
  my($lon,  $lat)	= lonLat($extent);
  my($colTo,$rowTo)	= flow_to($flowDir);
  my $up_area		= byte($flowDir)->upstrAccumAll(cell_area($lon,$lat,$extent));
  my($area, $cell_pdl)	= cell_table($up_area,$flowDir,$colTo,$rowTo);
  my @area		= $area->list;
  my $cell_table	= repack_data([$cell_pdl->list], [scalar(@area),5]);
  map push(@{$$cell_table[$_]}, Geo::GDAL::ApplyGeoTransform($$extent{gTransform},
	     $$cell_table[$_][3]+0.5, $$cell_table[$_][4]+0.5), $area[$_]), 0..$#$cell_table;
  push @time, Benchmark->new;
  printf "   %s\n",timestr(timediff($time[-1], $time[-2])) if $verbose;

		### Save data to files
  write_cellTable($file_ctNet,  ['ToCell','ToCellX','ToCellY','X','Y','Lon','Lat','SubbasinArea'],$cell_table);
  write_gridascii($file_upArNet,$up_area, $extent, {FORMAT => '%.1f'});
  write_gridascii($file_bID,    $bsnIDNet,$extNet, {FORMAT => '%d'});
  writeflex(	  $file_ctbNet, $up_area, $area, $cell_pdl);
}

#######################################################################
#######################################################################
						# Report Total Time
printf "\nTime used to build Cell Table - %d hours, %d minutes, and %d seconds\n\n",
	time_used($time_start,time()) if $verbose;

close OLDERR;
exit;

#######################################################################
######################  Functions  ####################################

sub trim

{
  my($net,$ext,$box,$tr) = @_;

  if ($tr == 1) {
    $$ext{ncols}     = $$box[2] - $$box[0] + 1;
    $$ext{nrows}     = $$box[1] - $$box[3] + 1;;
    $$ext{xllcorner}+= $$box[0]*$$ext{cellsize};
    $$ext{yllcorner}+= ($net->dim(1)-1-$$box[1])*$$ext{cellsize};

    $net	= $net($$box[0]:$$box[2],$$box[3]:$$box[1]);
  }
  if ($tr == 2) {
    my @trim	= (0,0,0,0);		### (top, bottom, left, right)
    my @check	= ($net->sumover, $net(:,-1:0)->sumover,
		   $net->transpose->sumover, $net(-1:0,:)->transpose->sumover);

    foreach my $i (0..3) { $trim[$i]++ while ($check[$i]->($trim[$i])->isbad) }
    $$ext{ncols}    -= $trim[2]+$trim[3];
    $$ext{nrows}    -= $trim[0]+$trim[1];
    $$ext{xllcorner}+= $trim[2]*$$ext{cellsize};
    $$ext{yllcorner}+= $trim[1]*$$ext{cellsize};

    $net	= $net($trim[2]:-$trim[3]-1,$trim[0]:-$trim[1]-1);
  }

  return $net, $ext;
}

#######################################################################

sub internal_basins

{
  my ($basinID, $flowDir) = @_;
  my  $basinInt	= zeroes($basinID->dims);

  my @mouthPnt	= whichND($flowDir == 0)->list;
  for (my $i=0; $i<$#mouthPnt; $i=$i+2) {
    if ($flowDir->range([$mouthPnt[$i]-1,$mouthPnt[$i+1]-1],3,'truncate')->ngood == 9) {
      my $id	= $basinID->at($mouthPnt[$i],$mouthPnt[$i+1]);
      $basinInt += $id * ($basinID == $id);
    }
  }

  return $basinInt->setvaltobad(0);
}

#######################################################################

sub write_cellTable

{
  my ($file, $header, $table) = @_;

  open (FILE,">$file") or die "Couldn't open $file, $!";
    print FILE join("\t",@$header),"\n";
    foreach my $row (@$table) {
      print FILE join("\t",@$row),"\n";
    }
  close FILE;
}

#######################################################################

# sub check_circularity_pdl
#
# {
#   my ($flowDir, $colTo, $rowTo, $gT) = @_;
#
#   my @dims = $flowDir->dims;
#   my $mask = zeroes(byte,@dims);
#
#   for (my $row=0; $row<$dims[1]; $row++) {
#     CC_LOOP:
#     for (my $col=0; $col<$dims[0]; $col++) {
#       next CC_LOOP if $flowDir($col,$row)->isbad || $flowDir($col,$row)==0 || $mask($col,$row)==1;
#       $mask($col,$row) .= 1;
#
#       my $col_to	= $colTo->at($col,$row);
#       my $row_to	= $rowTo->at($col,$row);
#
# # printf "Start= %d\t%d\t%d\t%d\t%d\t%d\n",$col,$row,$col_to,$row_to,$flowDir->at($col,$row),$flowDir->at($col_to,$row_to);
#       my $route_to	= $col_to.'_'.$row_to;
#       my %route		= ($col.'_'.$row => [$col,$row]);
#       my @flow		= ([$col,   $row,   $flowDir->at($col,   $row)],
# 			   [$col_to,$row_to,$flowDir->at($col_to,$row_to)]);
#       do {
# 	if (defined $route{$route_to}) {		### Circularity check
# 	  shift @flow while ($flow[0][0].'_'.$flow[0][1] ne $route_to);
# 	  my $route	= join("\n\t",map(join("\t",@$_),@flow));
# 	  my ($lon,$lat)= Geo::GDAL::ApplyGeoTransform($gT,$col_to+0.5,$row_to+0.5);
#
# 	  die <<END;
# 	Network circularity is detected at-
# 	(col,row) = ($col_to, $row_to)
# 	(lon,lat) = ($lon, $lat)
#
# 	Circular route-
# 	Col	Row	Dir
# 	$route
#
# 	Aborting...\n
# END
# 	}
# 	$route{$route_to}	= [$col_to,$row_to];
# 	my ($col_fr,$row_fr)	= ($col_to,$row_to);
# 	next CC_LOOP if $mask($col_fr,$row_fr)==1;
# 	$mask($col_fr,$row_fr) .= 1;
# 	push @flow, [$col_fr,$row_fr,$flowDir->at($col_fr,$row_fr)];
#
# 	$col_to		= $colTo->at($col_fr,$row_fr);
# 	$row_to		= $rowTo->at($col_fr,$row_fr);
# 	$route_to	= $col_to.'_'.$row_to;
# 	push @flow, [$col_to, $row_to, $flowDir->at($col_to,$row_to)];
# # print "Contd= ",$col_fr,"\t",$row_fr,"\t",$col_to,"\t",$row_to,"\t",$flowDir->at($col_fr,$row_fr),"\t",$flowDir->at($col_to,$row_to),"\n";
#       } while ( $flowDir($col_to,$row_to)>0);
#     }
#   }
# }

#######################################################################

sub script_dir
{
  $0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;
  return $1 || "./";
}

#######################################################################

sub usage

{
  my $app_name = basename($0);
  print <<EOF;

Usage:
	$app_name [-h] [-v] [-sub METHOD] [-ib] [-b BASIN_ID_FILE] [-p PROJ] NETWORK_FILE [CELL_TABLE_FILE]

This code generates CELL_TABLE_FILE (must be *.dat) from NETWORK_FILE (direction of flow).

Options:

h	Display this help.
v	Verbose mode.
sub	Make network subset by given method: <0- none; 0- set of basin IDs; 1- extents; 2- polygon mask.
ib	In addition, add internal basins to output.
b	Basin ID file. If omitted, the basin IDs will be ordered by watershed area.
p	Projection (default is epsg:4326).

EOF
  exit;
}

#######################################################################
###################  PDL::PP Functions  ###############################

__DATA__
__Pdlpp__

#######################################################################

pp_addhdr('
  #include <unistd.h>       /* we need defs of XXXX */
  #include <stdio.h>

  static void upstr(PDL_Byte * flowDir, PDL_Long * stack, PDL_Long n_size, PDL_Long m_size, PDL_Long N, PDL_Long M)
  {
    long from[2][8] = { {1,1,0,-1,-1,-1,0,1} , {0,1,1,1,0,-1,-1,-1} };
    long dir, ind, xx, yy;
    long NN = N;	long count = 1;
    long MM = M;	long pos   = 0;
    stack[0] = N + M*n_size;

    while (pos < count) {
      MM = stack[pos] / n_size;
      NN = stack[pos] - n_size*MM;
      pos++;
      for (dir=0; dir<8; dir++) {
        xx  = NN - from[0][dir];
        yy  = MM - from[1][dir];
        if (xx<0 || yy<0 || xx==n_size || yy==m_size) break;	// This never should happen
        ind = xx + yy*n_size;
        if (flowDir[ind] == (0x01<<dir)) stack[count++] = ind;
      }
    }
  }
');

#######################################################################

pp_def('upstreamMask', HandleBad => 1,
  Pars => 'byte flowDir(n,m);
    int N(); int M();
    byte [o] mask(n,m);',
  Code => '
    int ind;
    int n_size = $SIZE(n);    int NN = $N();
    int m_size = $SIZE(m);    int MM = $M();
    int *myStack;	myStack = malloc(n_size*m_size*sizeof *myStack);

    loop(n,m) %{ $mask() = 0; %}  	//	Initialization of the output arrays
    for (ind=0; ind < n_size*m_size; ind++) {
      myStack[ind] = -1;
    }

    upstr($P(flowDir),myStack,n_size,m_size,NN,MM);
    ind = 0;
    while (myStack[ind] != -1 && ind < n_size*m_size) {
      MM = myStack[ind] / n_size;
      NN = myStack[ind] - n_size*MM;
      ind++;
      $mask(n=>NN,m=>MM) = 1;
    }
');

#######################################################################

pp_def('upstrAccumAll', HandleBad => 1,
  Pars => 'byte flowDir(n,m); double data(n,m);
    double [o] upstrAccData(n,m);',
  Code => '
    int i, j, ind;
    int n_size = $SIZE(n);
    int m_size = $SIZE(m);
    double *myData;	myData  = malloc(n_size*m_size*sizeof *myData);
    int    *myStack;	myStack = malloc(n_size*m_size*sizeof *myStack);

	//	Initialization of the output arrays
    loop(n,m) %{ $upstrAccData() = 0; %}
    for (j=0; j<m_size; j++) {
      for (i=0; i<n_size; i++) {
	ind = i + j*n_size;
	myData[ind]  = ( $ISBAD($data(n=>i,m=>j)) ) ? 0 : $data(n=>i,m=>j);
	myStack[ind] = -1;
      }
    }

    for (j=0; j<m_size; j++) {
      for (i=0; i<n_size; i++) {
	if ( $ISBAD($flowDir(n=>i,m=>j)) ) {
	  $SETBAD($upstrAccData(n=>i,m=>j));
	}
	else {
	  upstr($P(flowDir),myStack,n_size,m_size,i,j);
	  ind = 0;
	  while (myStack[ind] != -1 && ind < n_size*m_size) {
	    $upstrAccData(n=>i,m=>j) += myData[myStack[ind]];
	    myStack[ind++] = -1;				// Reset myStack
	  }
	}
      }
    }
');

#######################################################################

pp_def('basinIDs', HandleBad => 1,
  Pars => 'byte flowDir(n,m);
    int table(k,l);
    int [o] basinID(n,m);',
  Code => '
    int ID     = 0;
    int n_size = $SIZE(n);
    int m_size = $SIZE(m);
    int l_size = $SIZE(l);
    int ind, i, j, i_cell;
    int *myStack;	myStack = malloc(n_size*m_size*sizeof *myStack);

	//	Initialization of the output arrays
    loop(n,m) %{ $basinID() = 0; %}
    for (j=0; j<m_size; j++) {
      for (i=0; i<n_size; i++) {
	if ( $ISBAD($flowDir(n=>i,m=>j)) ) $SETBAD($basinID(n=>i,m=>j));
	ind = i + j*n_size;
	myStack[ind] = -1;
      }
    }

	//	Building basin IDs
    for(i_cell=l_size-1; i_cell>=0; i_cell--) {		// Reverse order
      if ($table(k=>0,l=>i_cell) != 0) continue;
      i = $table(k=>3,l=>i_cell);
      j = $table(k=>4,l=>i_cell);

      upstr($P(flowDir),myStack,n_size,m_size,i,j);
      ind = 0;	ID++;
      while (myStack[ind] != -1 && ind < n_size*m_size) {
	j = myStack[ind] / n_size;
	i = myStack[ind] - n_size*j;
	myStack[ind++] = -1;				// Reset myStack
	$basinID(n=>i,m=>j) = ID;
      }
    }
');

#######################################################################

pp_def('checkCircularity', HandleBad => 1,
  Pars => 'byte flowDir(n,m);
    int [o] indexX();	int [o] indexY();',
  Code => '
    int n_size = $SIZE(n);
    int m_size = $SIZE(m);

	//	Initialization of "downstream" variables

    int i, j, xx, yy, dir, pix;
    int from[2][255];
    int From[2][9]	= { {0,1,1,0,-1,-1,-1,0,1} , {0,0,1,1,1,0,-1,-1,-1} };
    int dInd[9]		= {0,1,2,4,8,16,32,64,128};
    for (i=0; i<9; i++) {
      from[0][dInd[i]]	= From[0][i];
      from[1][dInd[i]]	= From[1][i];
    }
    $indexX() = -1;	$indexY() = -1;

	//	Initialization of myMask array

    int *myMask[n_size];		// Dynamically allocate memory
    for (i=0; i<n_size; i++)
      myMask[i] = (int *)malloc(m_size * sizeof(int));

    for (j=0; j<m_size; j++)		// Initialize with values
      for (i=0; i<n_size; i++)
	myMask[i][j] = -1;

	//	Calculations

    for (j=0; j<m_size; j++) {
      for (i=0; i<n_size; i++) {
	if ( $ISBAD($flowDir(n=>i,m=>j)) ) continue;

	pix = i + j * n_size;
	xx  = i;
	yy  = j;

	do {
	  if (myMask[xx][yy] > 0) break;		// Already checked pixel and flow downstream

	  myMask[xx][yy] = pix;			// Mark as checked pixel
	  dir = $flowDir(n=>xx,m=>yy);
	  xx += from[0][dir];
	  yy += from[1][dir];

	  if (dir && myMask[xx][yy] == pix) {	// Circularity found!!!
	    $indexX() = xx;	$indexY() = yy;
	    j = m_size;	i = n_size;	break;
	  }
	} while (myMask[xx][yy] == -1);
      }
    }
');

#######################################################################

pp_done();
