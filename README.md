# Supplementary scripts for meteorological equipment from RPG

The main aim of this repository is to make selected modules developed at RPG (Radiometer Physics GmbH) available for scietific (meteorological) community.

## License

Copyright (c) Radiometer Physics GmbH. All rights reserved.

The software is licensed under [GNU Lesser General Public License v2.1](./LICENSE).

## Acknowledgement

If you find the provided modules useful for you study, we kindly ask to acknowledge in publications our efforts to provide *free* and *open-source* solutions for the scientific community. These acknowledgements will help us to fulfil requirements for a publication in 
the Journal of Open Source Software.

## Found a bug?

Please contact using email: <remotesensing-service@radiometer-physics.de>

## Interested in contribution to the develoment of modules? Have an idea how to further improve functionality? 

Please send us your ideas and suggestions to <remotesensing-service@radiometer-physics.de>

## Included modules

### RadarControl.py

This module contains a Python API class that provides a platform independent interface to control RPG FMCW cloud radars. The class can be used for a development of adaptive observation strategies. The class provides an access to the latest measured sample (both integrated and spectral measured quantities). Based on user anaysis of the last sample, scanning and measurement settings can be switched in near-real time. The class can also be used to access housekeeping data and to monitor the radar status. In addition to the functions which can be used in the users scritps, the module also provides command line interface to control RPG FMCW radars directly from a command promt.

Usage of the module from a command prompt:

The Python module provides an assistance in using available commands. This way users do not need to browse through the code to understand how to use the class. In order to use the module please open the command promt and change the working directory to the one containing the python module. In order to get a list of available functions execute the following command in the command prompt (syntaxis is for Windows, in Linux and MAC OS it may differ):
```
> python RadarControl.py
```

In order to get a syntax template for a command (COMMAND) please execute the following:
```
> python RadarControl.py COMMAND "help"
```

Currently implemented functionality for the comand prompt use is listed below:

`get_radar_status` provides basic information about the current radar activity.
`start_radar_measurements` starts radar measurement using a MDF file on the host PC.
`start_radar_measurements_local_mdf` starts radar measurements using a MDF located on the user PC (MDF is not copied to the host).
`terminate_radar_measurements` stops currently running radar measurements.
`get_mdf_list` provides a list of MDF files available on the host PC in the default folder.
`get_radar_id` provides information about the radar.
`install_local_mdf` copies an MDF from the user PC to the host PC. The MDF is not started.