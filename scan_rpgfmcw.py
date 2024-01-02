# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 13:44:33 2023

@author: Toshiba
"""
import os
import time
import sys
import glob
import datetime
import numpy as np

from RPGtools_RemoteSensingInstruments.RadarControl import (Scan,
                                                            Client,
                                                            MeasDefFile,
                                                            MeasBatchFile,
                                                            install_local_mdf,
                                                            get_radar_status,
                                                            start_radar_measurements_local_mdf,
                                                            start_radar_measurements,
                                                            install_local_mdf)

WORKDIR = os.path.abspath('.')
WORKDIR = WORKDIR if WORKDIR.endswith(os.sep) else WORKDIR + os.sep

##################################
### RADAR SPECIFIC SETTINGS ######
##################################
IP = '192.168.0.2'
PORT = 7000
PW = ''
CONFIG = [IP, PORT, PW]

##################################
### SCAN SPECIFIC SETTINGS ###
##################################
SCANSPEED = 1
FASTSPEED = 5

##################################
### CHIRP SPECIFIC SETTINGS ######
##################################
CHIRPPRG = 7  # RMBLCHIRP, the inface number starts at 1, the real list at 0


##################################
### LOCATION SPECIFIC SETTINGS ###
##################################
NORTHOFFSET = 22.3


##################################
### CLOUDLAB SPECIFIC SETTINGS ###
##################################
# default duration is in s
DURATION = 20 * 60
# EXTRAWAIT = 10 # an extra number of seconds to wait for the scan to be done
# at the end of a scan, wait this amount in seconds before sending a new scan cmd
AFTERSCANWAIT = 5

##################################
### RPG SPECIFIC SETTINGS ###
##################################
t0 = datetime.datetime(2001, 1, 1)


def report(client, duration, sleeptime=1):
    for i in np.arange((duration/sleeptime) + 1):
        _sample = client.get_last_sample()
        dts = (datetime.timedelta(seconds=_sample.samp_t),
               datetime.timedelta(seconds=_sample.samp_ms/10**6))
        st = t0 + dts[0] + dts[1]
        print('Time of sample:', st)
        print('Current elevation:',
              _sample.elv, _sample.inc_el, _sample.inc_elax)
        print('Current azimuth:',
              _sample.azm, _sample.inc_elax, )
        time.sleep(sleeptime)


def ensure_termination(client,
                       timeout=20,
                       retrytime=0.5,
                       quiet=True):
    res = -1
    cnt = 0
    while res != 1:
        if not quiet:
            print(f'Trying to terminate measurements (try {cnt+1})')
        res = client.terminate_radar_measurements()

        if res in [3, 4]:
            print('Zero calibration cannot be determinated, wait longer...')
            try:
                time.sleep(1)
            except KeyboardInterrupt:
                print('Waiting cancelled')
                return
        elif res in [5]:
            print('Transmitter calibration cannot be determinated, wait longer...')
            time.sleep(10)

        time.sleep(retrytime)
        cnt += 1
        _status = client.get_radar_status()
        # print('Current status after termination command:', _status.__dict__)
        if cnt * retrytime >= timeout:
            print(f'Measurement could not be terminated in {
                  cnt*retrytime} seconds, use GUI')
            return

    # to ensure there is enough time for the radar to react
    time.sleep(2)


def ensure_start(client,
                 file,
                 timeout=30,
                 retrytime=2,
                 quiet=True):

    res = -1
    cnt = 0
    client.terminate_radar_measurements()
    while res != 1:
        if not quiet:
            print(f'Trying to start measurements (try {cnt+1})')

        # if the file exists we assume it is a local file that we do once
        # else we assume it is on the radar in the default MDF/MBF directory
        if os.path.exists(file):
            # assume its local and we need to send it
            res = client.start_radar_measurements_local_mdf(file)
        else:
            # assume its on the radar
            res = client.start_radar_measurements(file)

        # if res == 2:
        #     ensure_termination(client, quiet=quiet)

        time.sleep(retrytime)
        cnt += 1
        _status = client.get_radar_status()
        if 'mdf_name' in _status.__dict__:
            curmdf = _status.mdf_name

            if isinstance(curmdf, list):
                curmdf = curmdf[0]

            if file.lower().endswith(curmdf.lower()):
                print('Radar reports matching MDF')
                return

        # print('Current status after termination command:', _status.__dict__)
        if cnt * retrytime >= timeout:
            print(f'Measurement could not be started in {
                  cnt} seconds, use GUI')
            return

    # to ensure there is enough time for the radar to react
    time.sleep(2)


def scan_rhi(elevation_init,
             elevation_end,
             azimuth,
             **kwargs):
    scan_generic(elevation_init, elevation_end, azimuth, azimuth,
                 once=True,
                 **kwargs)


def scan_ppi(elevation,
             **kwargs):
    azimuth_init, azimuth_end = 0 + NORTHOFFSET, 359.99 + NORTHOFFSET
    scan_generic(elevation, elevation, azimuth_init, azimuth_end,
                 once=True,
                 **kwargs)


def scan_elevation(elevation_init,
                   elevation_end,
                   azimuth,
                   **kwargs):

    scan_generic(elevation_init, elevation_end, azimuth, azimuth, **kwargs)


def scan_azimuth(azimuth_init,
                 azimuth_end,
                 elevation,
                 **kwargs):
    scan_generic(elevation, elevation, azimuth_init, azimuth_end, **kwargs)


def scan_generic(elevation_init,
                 elevation_end,
                 azimuth_init,
                 azimuth_end,
                 # at which speed to scan
                 scanspeed=SCANSPEED,
                 # at which speed to move
                 fastspeed=FASTSPEED,
                 # this is the LOWER scantime, as it will be updated
                 # to match an even number of scans with the given speed/angle
                 duration=DURATION,
                 # the calibration interval (abs. cal.). be aware that once
                 # this is running the scan cannot be aborted
                 calibration_interval=1,
                 reporting=True,
                 reportinterval=5,
                 once=False,
                 quiet=True,
                 dryrun=True):
    try:
        c = Client(*CONFIG, SuppressOutput=quiet)
    except:
        print('Error connecting to data server, returning ...')
        print('Is the CLIENT running (on this computer).',
              'We need the data server to forward commands ...')
        return

#
    if (azimuth_init == azimuth_end or elevation_init == elevation_end):
        if azimuth_init == azimuth_end:
            # rhi like scan, movement in azi is fast
            movementtime = abs(90 - elevation_init) / scanspeed
            movementtime += abs(azimuth_init-NORTHOFFSET) / fastspeed
        elif elevation_init == elevation_end:
            # sector scan/ppi like scan, movement in elv is fast
            movementtime = abs(90 - elevation_init) / fastspeed
            movementtime += abs(azimuth_init - NORTHOFFSET) / scanspeed

        movementtime = int(np.ceil(movementtime))
        print(movementtime, '**')
    else:
        print('Scanning in both azimuth and elevation is not supported by',
              'this script. exiting....')
        return

    if azimuth_init == azimuth_end:
        onescanduration = (abs(elevation_end - elevation_init)/scanspeed)
        onescanduration = int(onescanduration)
    elif elevation_init == elevation_end:
        onescanduration = (abs(azimuth_end - azimuth_init)/scanspeed)
        onescanduration = int(onescanduration)
    else:
        print('Scanning in both azimuth and elevation is not supported by',
              'this script. exiting....')
        return

    if once:
        duration = onescanduration + movementtime
        nscans = 1
    else:
        nscans = duration / onescanduration
        nscans = int(np.ceil(nscans))
        # these have to be symmetrical, so always needs to be an even number
        if nscans % 2 != 0:
            print('Adding another scancycle to return radar to zenith')
            nscans += 1

        duration = int(nscans * onescanduration) + movementtime

    print(f'The scanrange of {elevation_init}° to {elevation_end}° with',
          f'a speed of {scanspeed} results in',
          f' {nscans} scans for {duration} seconds',
          '(This was rounded up to match the next higher duration for n scans)')

    if azimuth_init == azimuth_end:
        mdffilename = 'ELEVATION_SCAN.MDF' if not once else 'RHI_SCAN.MDF'
        # elevation scan
        # the first scan going down
        SCAN_FORTH = Scan(elv=elevation_init,
                          azm=((azimuth_init-NORTHOFFSET)+360) % 360,
                          elv_target=elevation_end,
                          azm_target=((azimuth_init-NORTHOFFSET)+360) % 360,
                          elv_spd=scanspeed,
                          azm_spd=fastspeed,
                          )
        # the second scan (e.g. going back up)
        SCAN_BACK = Scan(elv=elevation_end,
                         azm=((azimuth_init-NORTHOFFSET)+360) % 360,
                         elv_target=elevation_init,
                         azm_target=((azimuth_init-NORTHOFFSET)+360) % 360,
                         elv_spd=scanspeed,
                         azm_spd=fastspeed,
                         )
    elif elevation_init == elevation_end:
        # maybe a bit of a misnomer but follows the CLOUDLAB nomenclature
        # from the miras
        mdffilename = 'SECTOR_SCAN.MDF' if not once else 'PPI_SCAN.MDF'
        # azimuth/sector scan
        # the first scan going down
        SCAN_FORTH = Scan(elv=elevation_init,
                          azm=((azimuth_init-NORTHOFFSET)+360) % 360,
                          elv_target=elevation_init,
                          azm_target=((azimuth_end-NORTHOFFSET)+360) % 360,
                          elv_spd=fastspeed,
                          azm_spd=scanspeed,
                          )
        # the second scan (e.g. going back up)
        SCAN_BACK = Scan(elv=elevation_init,
                         azm=((azimuth_end-NORTHOFFSET)+360) % 360,
                         elv_target=elevation_init,
                         azm_target=((azimuth_init-NORTHOFFSET)+360) % 360,
                         elv_spd=fastspeed,
                         azm_spd=scanspeed,
                         )

    SCANS = [SCAN_FORTH, SCAN_BACK]
    # once means a RHI or PPI scan
    if once:
        SCANS = SCANS[:1]

    m = MeasDefFile()

    if once:
        frames = [[0, 0, 1]]
    else:
        frames = [[0, 1, int(np.ceil(nscans/2))]]

    if dryrun:
        m.create(WORKDIR + mdffilename,
                 CHIRPPRG,
                 SCANS,
                 frames=frames,
                 duration=duration,
                 filelen=onescanduration,
                 cal_int=calibration_interval,
                 )
        m.read(WORKDIR + mdffilename)
        print(f'Made {WORKDIR+mdffilename}:')
        m.output()
        # os.remove(WORKDIR + f)
    else:
        radar_id = c.get_radar_id()
        radar_status = c.get_radar_status()
        if 'mdf_name' in radar_status.__dict__:
            oldmdf = radar_status.mdf_name
        else:
            oldmdf = None

        if isinstance(oldmdf, list):
            oldmdf = oldmdf[0]

        print(f'Radar is currently running {oldmdf}\n')
        try:
            print(f'A Running {WORKDIR+mdffilename}:')
            # m.create(WORKDIR + f, CHIRPPRG, SCAN_FORTH, duration=duration)
            m.create(WORKDIR + mdffilename,
                     CHIRPPRG,
                     SCANS,
                     frames=frames,
                     duration=duration,
                     filelen=onescanduration,
                     cal_int=calibration_interval,
                     )
            m.read(WORKDIR + mdffilename)
            ensure_start(c, WORKDIR+mdffilename)
            time.sleep(max([5, movementtime]))
            time.sleep(10)
            if reporting:
                report(c, duration + 1, sleeptime=reportinterval)
            else:
                time.sleep(duration + 1)

            time.sleep(AFTERSCANWAIT)
        except KeyboardInterrupt:
            print('Stopping scanning operation...')
        finally:
            c.terminate_radar_measurements()

        print('Scan finished (see above if successful).')

        if oldmdf is not None:
            print(f'Installing previous MDF: {oldmdf}')
            # time.sleep(15)
            # ensure_termination(c)
            time.sleep(10)
            # c.terminate_radar_measurements()
            ensure_start(c, oldmdf)

        # ensure some wait before killing the client when not in dryrun
        time.sleep(5)
    return c
    # del c


if __name__ == '__main__':
    pass
    # a half rhi at positioner 0° to avoid unneccessary movement for testing
    client = scan_elevation(90, 60,
                            NORTHOFFSET,
                            duration=20,
                            # dryrun=False,
                            dryrun=True,
                            quiet=True)

    client = scan_ppi(85,
                       dryrun=False,
                      # dryrun=True,
                      quiet=False)

    client = scan_rhi(10, 170, NORTHOFFSET,
                      # dryrun=False,
                       dryrun=True,
                      quiet=True)

    # make some MBF file from all MDFs you have.
    # mbf = MeasBatchFile()
    # mdflist = [WORKDIR + i for i in os.listdir(WORKDIR) if 'MDF' in i]
    # mbf.create(WORKDIR+'test.mbf', mdflist, repetitions=3)
