# Supplementary scripts for meteorological equipment from RPG

The main aim of this repository is to make selected modules developed at RPG (Radiometer Physics GmbH) available for the scientific (meteorological) community.

## License

Copyright (c) Radiometer Physics GmbH. All rights reserved.

The software is licensed under [GNU Lesser General Public License v2.1](./LICENSE).

## Acknowledgement

If you find the provided modules useful for you study, we kindly ask to acknowledge in publications our efforts to provide *free* and *open-source* solutions for the scientific community. These acknowledgements will help us to fulfil requirements for a publication in 
the Journal of Open Source Software.

## Found a bug?

Please contact us via email: <remotesensing-service@radiometer-physics.de>

## Interested in contributing to the develoment of modules? Have an idea how to further improve functionality? 

Please send us your ideas and suggestions to <remotesensing-service@radiometer-physics.de>

## Included modules

### RadarControl.py

This module contains a Python API class that provides a platform independent interface to control RPG FMCW cloud radars. The class can be used for a development of adaptive observation strategies. The class provides an access to the latest measured sample (both integrated and spectral measured quantities). Based on user anaysis of the last sample, scanning and measurement settings can be switched in near-real time. The class can also be used to access housekeeping data and to monitor the radar status. In addition to the functions which can be used in the users scripts, the module also provides command line interface to control RPG FMCW radars directly from a command prompt.

Usage of the module from a command prompt:

The Python module provides an assistance in using available commands. This way users do not need to browse through the code to understand how to use the class. In order to use the module, please open the command promt and change the working directory to the one containing the python module. In order to get a list of available functions, execute the following command in the command prompt (syntaxis is for Windows, in Linux and MAC OS it may differ):
```
> python RadarControl.py
```

In order to get a syntax template for a command (COMMAND) please execute the following:
```
> python RadarControl.py COMMAND "help"
```

Currently implemented functionality for the command prompt use is listed below:

`get_radar_status` provides basic information about the current radar activity.

`start_radar_measurements` starts radar measurement using a MDF file on the host PC.

`start_radar_measurements_local_mdf` starts radar measurements using a MDF located on the user PC (MDF is not copied to the host).

`terminate_radar_measurements` stops currently running radar measurements.

`get_mdf_list` provides a list of MDF files available on the host PC in the default folder.

`get_radar_id` provides information about the radar.

`install_local_mdf` copies an MDF from the user PC to the host PC. The MDF is not started.

### scan_rpgfmcw.py

This module contains utilities/examples to conduct scans. The following scanpatterns are implements/available:
- PPI scan with a mandatory elevation input
- RHI scan with mandatory azimuth, and elevation (init and end)
- repeated RHI scans similar to the above RHI scan with a duration keyword. Generally, a forth and a back RHI are conducted, i.e. the number of scans is even
- repeated partial PPI (sector/azimuth) scans with a given elevation and azimuth (init and end) for a selectable duration

The former two are standard scans and should not require additional information. The latter two are meant to be repeatable (and reduced for the sector scans) RHI/PPI scans to better track changes over time. Generally, the `scan_generic` function is the central function and all other functions are shallow wrappers around it as the logic is similar for all scans it is located in `scan_generic` rather than the single scans. Each scan returns the original Client from RadarControl.py to allow changes afterwards. 

**At the start of the file several configuration parameters are available that should be checked and adjusted, especially the scan speeds and the north angle**

The working principle is to create a MDF files with a set number of frame repetitions (number of scans that fit in duration / 2) which is then sent to the data server (which communicates with the radar). After a scan has been conducted, the previous MDF file is reinstalled (no support for MBF reinstallation yet). Dryrun creates MDF files without starting/sending them, e.g. for error checking. Dryrun is set to True per default to allow checking of the approach before.

#### Example
```
# set your own for the specific location, this should be set 
NORTHOFFSET = 0
# Does a rhi scan of 30° (for illustration) for 20 seconds
# (this gets rounded up to 30 seconds to fit the anglerange (30° movement) with
# a default scanspeed of 1°/s)
client = scan_elevation(90, 60,
                        NORTHOFFSET,
                        duration=20,
                        dryrun=True,
                        quiet=True)

# PPI scan, simply choose the elevation at which to conduct
client = scan_ppi(85,
                  dryrun=True,
                  quiet=False)

# RHI scan from 10 to 170° elevation at the NORTH offset
client = scan_rhi(10, 170, NORTHOFFSET,
                  dryrun=True,
                  quiet=True)
```

#### Issues
Currently, the termination/ensuring start of the scan sometimes faces issues, when 
- the original MDF has been restarted and has a long zero calibration value a repetition of the scan won't work
- the reporting of the values during the scanning sometimes indicates weird time stamps, mainly when for some reason the same MDF gets restarted.

