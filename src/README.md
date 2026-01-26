Workflow of `src.py`
```FOR each tropical cyclone
    Load track and define simulation time window
    Load storm tide, reanalysis, wind, and rainfall data
    Align all datasets temporally
    Calculate tide timing offsets
    Merge storm tide with shifted reanalysis tides
    Run or load SFINCS simulations
    Extract maximum hazard metrics
    Compute flood depth
    Attribute flooding to coastal, runoff, or compound processes
END```
