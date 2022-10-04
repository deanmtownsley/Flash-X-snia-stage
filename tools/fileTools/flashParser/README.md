### Parser for Flash-X log files to collect performance data and statistics

Installation in development mode
```
python3 setup.py develop --user
```

Installation in produciton mode
```
python3 setup.py install
```

Usage works with with both `.log` and `.log.csv` files
```
import flashparser

log_dict = flashparser.LogDict("/path/to/logfile")

# Get number of calls for a timer
log_dict["Grid_updateRefinement"]["num calls"]

# Average time per proc
log_dict["Grid_updateRefinement"]["avg/proc"]

# Time per processor
log_dict["Grid_updateRefinement"]["proc/0000"]
log_dict["Grid_updateRefinement"]["proc/0001"]
.
.
log_dict["Grid_updateRefinement"]["proc/0010"]
.
.

# Get timer level
log_dict["Grid_updateRefinement"]["level"]
```
