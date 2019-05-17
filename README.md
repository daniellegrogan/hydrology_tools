# hydrology_tools
Tools for WBM work flow and processing hydrology data written by Alex Prusevich (alex.proussevitch@unh.edu)

Files included:

**spatial_aggregation.pl**
This script aggregates gridded values over areas based on an input mask and associated ID file.  Default aggregation is to average all values. Use -h flag to access help information and options.

**temporal_aggregation.pl**
This script aggregates gridded values by time, with options to average by month, year, or (daily, monthly, yearly) climatology. Default aggregation is average value over time.  Use -h flag to ccess help information and options.

**cell_table_pp.pl**
This script subsets river networks and generates files for:
- Network direction
- Basin IDs
- Upstream area
See cell_table_manual.init for a commented example of an input file and options to use this script.
