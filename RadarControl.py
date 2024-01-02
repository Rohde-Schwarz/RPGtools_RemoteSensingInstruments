# This module contains a Python API class that provides a platform independent interface
# to control RPG FMCW cloud radars. The class can be used for a development of adaptive
# observation strategies. The class provides an access to the latest measured sample
# (both integrated and spectral measured quantities). Based on user anaysis of the last
# sample, scanning and measurement settings can be switched in near-real time. The class
# can also be used to access housekeeping data and to monitor the radar status.
# In addition to the functions which can be used in the users scripts, the module also
# provides command line interface to control RPG FMCW radars directly from a command prompt.
#
# LICENSE
#
# Copyright (c) Radiometer Physics GmbH. All rights reserved.
#
# The software is licensed under GNU Lesser General Public License v2.1.
#
# ACKNOWLEDGEMENT
#
# If you find the provided modules useful for you study, we kindly ask
# to acknowledge in publications our efforts to provide *free* and
# *open-source* solutions for the scientific community. These acknowledgements
# will help us to fulfil requirements for a publication in the Journal
# of Open Source Software.
#
# Found a bug? Please contact us via email: <remotesensing-service@radiometer-physics.de>
#
# Interested in contributing to the development of modules? Have an idea how to further
# improve functionality? Please send us your ideas and suggestions to
# <remotesensing-service@radiometer-physics.de>

import socket
import sys
import os
import mmap
import struct
import numpy as np
import math

from array import array
from pathlib import Path


class ByteReader:

    def __init__(self, byte_input, start_idx=1):

        # index starts from 1 because the 0th byte
        # in the host response repeats the requested command
        # This repeated byte is checked
        # in the __check_first_byte_of_response fucntion
        # of the Client class (see below in this module)
        self.__i = start_idx
        self.__input = byte_input

    def __get_byte_num(self, value_type):

        # float, integer, unsigned integer
        if value_type in ('f', 'i', 'I'):
            return 4  # [bytes]

        # byte
        if value_type == 'b':
            return 1  # [byte]

        # short integer, unsigned short integer
        if value_type in ('h', 'H'):
            return 2  # [bytes]

        if value_type in ('q', 'Q'):
            return 8  # [bytes]

        raise Exception("Unknown variable type")

    def __read_single_value(self, value_type, endian='<'):

        V = self.__read_vector(value_type, 1, endian=endian)
        return V[0]

    def __read_vector(self, value_type, value_num, endian='<'):

        size = self.__get_byte_num(value_type)
        firstidx = self.__i
        lastidx = self.__i + value_num * size

        F = struct.unpack(value_type * value_num,
                          self.__input[firstidx:lastidx])

        self.__i += value_num * size
        return F

    def _read_float(self):

        try:
            return self.__read_single_value('f')
        except:
            raise Exception('ByteReader: Cannot read a float')

    def _read_unsigned_short(self):

        try:
            return self.__read_single_value('H')
        except:
            raise Exception('ByteReader: Cannot read a unsigned short int')

    def _read_signed_short(self):

        try:
            return self.__read_single_value('h', '<')
        except:
            raise Exception('ByteReader: Cannot read a signed short int ')

    def _read_unsigned_int(self):

        try:
            return self.__read_single_value('I')
        except:
            raise Exception('ByteReader: Cannot read a unsigned int')

    def _read_int(self):

        try:
            return self.__read_single_value('i')
        except:
            raise Exception('ByteReader: Cannot read an signed integer')

    def _read_int_big_endian(self):

        try:
            return self.__read_single_value('i', '>')
        except:
            raise Exception('ByteReader: Cannot read an integer')

    def _read_int_little_endian(self):

        try:
            return self.__read_single_value('i', '<')
        except:
            raise Exception('ByteReader: Cannot read an integer')

    def _read_long_long(self):

        try:
            return self.__read_single_value('q')
        except:
            raise Exception('ByteReader: Cannot read long long')

    def _read_byte(self):

        try:
            return self.__read_single_value('b')
        except:
            raise Exception('ByteReader: Cannot read a byte')

    def _read_string(self):

        try:
            for j in range(self.__i, len(self.__input)):
                if self.__input[j] == 0:
                    break

            s = self.__input[self.__i:j].decode('UTF-8')
            self.__i = j + 1
            return s
        except:
            raise Exception('ByteReader: Cannot read a string')

    def _read_null_separated_strings_to_list(self):

        try:
            a = self.__input[self.__i:].decode('UTF-8')
            b = a.split('\x00')
            return b
        except:
            raise Exception('ByteReader: Cannot read strings')

    def _read_float_vector(self, value_num):

        try:
            return self.__read_vector('f', value_num)
        except:
            raise Exception('ByteReader: Cannot read float vector')

    def _read_int_vector(self, value_num):

        try:
            return self.__read_vector('i', value_num)
        except:
            raise Exception('ByteReader: Cannot read int vector')

    def _read_byte_vector(self, value_num):

        try:
            return self.__read_vector('b', value_num)
        except:
            raise Exception('ByteReader: Cannot read byte vector')

    def _read_short_int_vector(self, value_num):

        try:
            return self.__read_vector('h', value_num)
        except:
            raise Exception('ByteReader: Cannot read short int vector')


class MDFList(ByteReader):

    def __init__(self, byte_input, SuppressOutput=False):

        self.__suppressout = SuppressOutput
        ByteReader.__init__(self, byte_input)
        self.__read_mdf_list()

    def __output(self):

        if self.__suppressout:
            return

        print('\nList of MDF files on the host PC: \n')

        for mdf in self.mdf_list:
            print(mdf)

    def __read_mdf_list(self):

        self.mdf_num = ByteReader._read_int(self)
        self.mdf_list = ByteReader._read_null_separated_strings_to_list(self)
        self.__output()


class Status(ByteReader):

    def __init__(self, byte_input, SuppressOutput=False):

        self.__suppressout = SuppressOutput
        ByteReader.__init__(self, byte_input)
        self.__read_radar_status()

    def __check_connection(self):

        self.connection = ByteReader._read_byte(self)

        if self.connection == 0:
            if self.__suppressout == False:
                print('The HOST is not connected to the radar')
            return False

        if self.connection != 1:
            if self.__suppressout == False:
                print('Invalid response: ', self.connection)
            return False

        if self.__suppressout == False:
            print('The HOST is connected to the radar')
        return True

    def __output_status(self):

        if self.__suppressout:
            return

        if self.connection == 0:
            print('The HOST is not connected to the radar')
            return

        if self.connection != 1:
            print('Invalid response')
            return

        if self.status == 1:
            print('The radar is in STANDBY mode')

        if self.status == 2:
            print('Measurement running')

        if self.status == 3:
            print('Zero calibration running')

        if self.status == 4:
            print('Absolute calibration running')

        if self.status == 5:
            print('Transmitter calibration running')

        if self.status == 2 or self.status == 3 or self.status == 5:

            print('Number of MDFs in current measurement: {}'.format(self.mdf_num))

            for i in range(self.mdf_num):

                print('MDF {}: '.format(i+1) + self.mdf_name[i])

            if self.mdf_num > 1:

                print('MBF: ' + self.mbf_name)
                print('Number of repetitions in batch: {}'.format(self.RepN))
                print('Current batch repetition index: {}'.format(self.CurRep))
                print('Current MDF index: {}'.format(self.CurMDF))

        if self.hierarchy == 0:
            print('Single radar')

        if self.hierarchy == 1:
            print('Dual radar, radar is in Master mode')

        if self.hierarchy == 2:
            print('Dual radar, radar is in Slave mode')

    def __read_radar_status(self):

        if self.__check_connection() == False:
            return

        self.status = ByteReader._read_byte(self)

        if self.status == 2 or self.status == 3 or self.status == 5:

            self.mdf_num = ByteReader._read_int(self)

            self.mdf_name = []

            for i in range(self.mdf_num):

                self.mdf_name.append(ByteReader._read_string(self))

            if self.mdf_num > 1:

                self.mbf_name = ByteReader._read_string(self)
                self.RepN = ByteReader._read_int(self)
                self.CurRep = ByteReader._read_int(self)
                self.CurMDF = ByteReader._read_int(self)

        self.hierarchy = ByteReader._read_byte(self)

        self.__output_status()


class RadarID(ByteReader):

    def __init__(self, byte_input, SuppressOutput=False):

        self.__suppressout = SuppressOutput
        ByteReader.__init__(self, byte_input)
        self.__read_radar_id()

    def __output(self):

        if self.__suppressout:
            return

        if self.connection == 0:
            print('The HOST is not connected to the radar')
            return

        if self.connection != 1:
            print('Invalid response: ', self.connection)
            return

        print('Software version: {0:.2f}'.format(self.swversion))
        print('Subversion: {}'.format(self.sbversion))

        if self.model == 0:
            print('Radar model: RPG-FMCW-94-SP')
        if self.model == 1:
            print('Radar model: RPG-FMCW-94-DP')
        if self.model == 2:
            print('Radar model: RPG-FMCW-35-SP')
        if self.model == 3:
            print('Radar model: RPG-FMCW-35-DP')

        print('Fabrication year: {}'.format(self.year))
        print('Fabrication number: {}'.format(self.number))
        print('Customer: ' + self.customer)
        print('License: {}'.format(self.license))

        if self.scanner == 0:
            print('No positioner found')
        if self.scanner == 1:
            print('Positioner found')

        if self.polmode == 0:
            print('Single pol. radar')
        if self.polmode == 1:
            print('Dual pol. radar: LDR mode')
        if self.polmode == 2:
            print('Dual pol. radar: STSR mode')

        print('IF range min. frequency [Hz]: {0:.1f}'.format(self.IFmin))
        print('IF range max. frequency [Hz]: {0:.1f}'.format(self.IFmax))
        print('Wavelength [m]: {0:.4f}'.format(self.wavelen))
        print('Antenna diameter [m]: {0:.3f}'.format(self.ant_d))
        print('Antenna separation [m]: {0:.3f}'.format(self.ant_sep))
        print('Antenna gain [linear]: {0:.3f}'.format(self.ant_g))
        print('Half-power-beam-width [deg]: {0:.3f}'.format(self.hpbw))
        print('Subreflector blockage [%]: {0:.3f}'.format(self.ant_block))
        print('Receiver gain V-pol. [linear]: {0:.3f}'.format(self.vrec_g))
        print('Receiver gain H-pol. [linear]: {0:.3f}'.format(self.hrec_g))

        if self.fft_win == 0:
            print('FFT window: Rectangular')
        if self.fft_win == 1:
            print('FFT window: Parzen')
        if self.fft_win == 2:
            print('FFT window: Blackman')
        if self.fft_win == 3:
            print('FFT window: Welch')
        if self.fft_win == 4:
            print('FFT window: Slepian2')
        if self.fft_win == 5:
            print('FFT window: Slepian3')

        if self.recovery == 0:
            print('Recovery after power failure: Disabled')
        if self.recovery == 1:
            print('Recovery after power failure: Enabled')

        if self.calibrat == 0:
            print('Absolute calibration: Not available')
        if self.calibrat == 1:
            print('Absolute calibration: Available')

        if self.pow_diag == 0:
            print('Power supply diagnostic: Not installed')
        if self.pow_diag == 1:
            print('Power supply diagnostic: Installed')

        print('Positioner azimuth offset [deg]: {0:.3f}'.format(self.scan_off))

        if self.interlan == 0:
            print('InterLAN status: Detection disabled')
        if self.interlan == 1:
            print('InterLAN status: Autodetection')

        print('Radar IP: ' + self.radar_ip)

    def __read_radar_id(self):

        if self.__check_connection() == False:
            return

        self.swversion = ByteReader._read_float(self)
        self.sbversion = ByteReader._read_int(self)
        self.model = ByteReader._read_int(self)
        self.year = ByteReader._read_int(self)
        self.number = ByteReader._read_int(self)
        self.customer = ByteReader._read_string(self)
        self.license = ByteReader._read_int(self)
        self.scanner = ByteReader._read_byte(self)
        self.polmode = ByteReader._read_byte(self)
        self.IFmin = ByteReader._read_float(self)
        self.IFmax = ByteReader._read_float(self)
        self.wavelen = ByteReader._read_float(self)
        self.ant_d = ByteReader._read_float(self)
        self.ant_sep = ByteReader._read_float(self)
        self.ant_g = ByteReader._read_float(self)
        self.hpbw = ByteReader._read_float(self)
        self.ant_block = ByteReader._read_float(self)
        self.vrec_g = ByteReader._read_float(self)
        self.hrec_g = ByteReader._read_float(self)
        self.fft_win = ByteReader._read_int(self)
        self.recovery = ByteReader._read_byte(self)
        self.calibrat = ByteReader._read_byte(self)
        self.pow_diag = ByteReader._read_byte(self)
        self.scan_off = ByteReader._read_float(self)
        self.interlan = ByteReader._read_byte(self)
        self.radar_ip = ByteReader._read_string(self)

        self.__output()

    def __check_connection(self):

        self.connection = ByteReader._read_byte(self)

        if self.connection == 0:
            if self.__suppressout == False:
                print('The HOST is not connected to the radar')
            return False

        if self.connection != 1:
            if self.__suppressout == False:
                print('Invalid response')
            return False

        if self.__suppressout == False:
            print('The HOST is connected to the radar')
        return True


class HouseKeeping(ByteReader):

    def __init__(self, byte_input, SuppressOutput=False):

        self.__suppressout = SuppressOutput
        ByteReader.__init__(self, byte_input)

        self.__read()

    def __read(self):

        self.gps_flag = ByteReader._read_byte(self)

        if self.gps_flag == 1:

            self.pos_stat = ByteReader._read_byte(self)
            self.time_stat = ByteReader._read_byte(self)

            if self.pos_stat == 1:

                self.longitude = ByteReader._read_float(self)
                self.latitude = ByteReader._read_float(self)
                self.gps_pos_time = ByteReader._read_int(self)

            if self.time_stat == 1:

                self.gps_sync_time = ByteReader._read_int(self)

        self.radar_time = ByteReader._read_int(self)
        self.met_found = ByteReader._read_byte(self)

        if self.met_found == 1:

            self.env_temp = ByteReader._read_float(self)
            self.pressure = ByteReader._read_float(self)
            self.rel_hum = ByteReader._read_float(self)
            self.wind_sp = ByteReader._read_float(self)
            self.wind_dir = ByteReader._read_float(self)

        self.scanner_found = ByteReader._read_byte(self)

        if self.met_found == 1:

            self.elev = ByteReader._read_float(self)
            self.azm = ByteReader._read_float(self)

        self.rec_temp = ByteReader._read_float(self)
        self.trans_temp = ByteReader._read_float(self)
        self.pc_temp = ByteReader._read_float(self)
        self.rain_stat = ByteReader._read_byte(self)

        self.heat_switch = ByteReader._read_int(self)
        self.blow_switch = ByteReader._read_int(self)
        self.rad_srive_cnt = ByteReader._read_int(self)

        self.free_mem = []
        self.tot_mem = []

        for i in range(self.rad_srive_cnt):

            self.free_mem.append(ByteReader._read_int(self))
            self.tot_mem.append(ByteReader._read_int(self))

        self.inc_el = ByteReader._read_float(self)
        self.inc_el_ax = ByteReader._read_float(self)
        self.meas_mode = ByteReader._read_byte(self)
        self.hirarchy = ByteReader._read_byte(self)

        if self.hirarchy == 1:

            self.sl_model_no = ByteReader._read_int(self)

        self.sl_error = ByteReader._read_byte(self)
        self.meas_run = ByteReader._read_byte(self)
        self.hdd_overflow = ByteReader._read_byte(self)
        self.alarm_code = ByteReader._read_byte(self)


class LastSample(ByteReader):

    def __init__(self, byte_input, SuppressOutput=False):

        self.__suppressout = SuppressOutput
        ByteReader.__init__(self, byte_input)

        self.__default_dr = 30   # [m]
        self.__defauls_dv = 0.04  # [m/s]

        self.__read_last_sample()

    def __check_running_measurements(self):

        self.meas_run = ByteReader._read_byte(self)

        if self.meas_run == 0:
            if self.__suppressout == False:
                print('Measurements are not running')
            return False

        if self.meas_run != 1:
            if self.__suppressout == False:
                print('Invalid response')
            return False

        if self.__suppressout == False:
            print('Measurements are running')
        return True

    def __read_last_sample(self):

        if self.__check_running_measurements() == False:
            return

        self.samp_idx = ByteReader._read_int(self)
        self.head_len = ByteReader._read_int(self)
        self.cgprog = ByteReader._read_int(self)
        self.model = ByteReader._read_int(self)
#         self.number    = ByteReader._read_int_big_endian(self)

        self.progname = ByteReader._read_string(self)
        self.customer = ByteReader._read_string(self)

        self.freq = ByteReader._read_float(self)
        self.ant_sep = ByteReader._read_float(self)
        self.ant_d = ByteReader._read_float(self)
        self.ant_g = ByteReader._read_float(self)
        self.hpbw = ByteReader._read_float(self)
        self.radar_c = ByteReader._read_float(self)
        self.pol_mode = ByteReader._read_byte(self)
        self.comp_ena = ByteReader._read_byte(self)
        self.antialia = ByteReader._read_byte(self)
        self.samp_dur = ByteReader._read_float(self)
        self.gps_lat = ByteReader._read_float(self)
        self.gps_lon = ByteReader._read_float(self)
        self.call_int = ByteReader._read_int(self)
        self.raltn = ByteReader._read_int(self)
        self.taltn = ByteReader._read_int(self)
        self.haltn = ByteReader._read_int(self)
        self.sequn = ByteReader._read_int(self)
        self.ralts = ByteReader._read_float_vector(self, self.raltn)
        self.talts = ByteReader._read_float_vector(self, self.taltn)
        self.halts = ByteReader._read_float_vector(self, self.haltn)
        self.nfft = ByteReader._read_int_vector(self, self.sequn)
        self.rng_off = ByteReader._read_int_vector(self, self.sequn)
        self.c_rep = ByteReader._read_int_vector(self, self.sequn)
        self.seqint = ByteReader._read_float_vector(self, self.sequn)
        self.dr = ByteReader._read_float_vector(self, self.sequn)
        self.max_vel = ByteReader._read_float_vector(self, self.sequn)
        self.center_f = ByteReader._read_float_vector(self, self.sequn)
        self.cal_par = ByteReader._read_int(self)
        self.samp_rate = ByteReader._read_int(self)
        self.max_range = ByteReader._read_int(self)

#         range_bin_num = math.ceil(self.max_range / self.__default_range_res)
#         vel_bin_num   = math.ceil(max(self.max_vel) / self.__defaulf_vel_res)

#         self.spectrum    = np.empty((vel_bin_num,range_bin_num))
#         self.spectrum[:] = np.nan

        self.sup_pow_lev = ByteReader._read_byte(self)
        self.spk_filter = ByteReader._read_byte(self)
        self.phase_corr = ByteReader._read_byte(self)
        self.rel_pow_cor = ByteReader._read_byte(self)
        self.fft_window = ByteReader._read_byte(self)

        self.fft_input_range = ByteReader._read_int(self)
        self.noise_filter = ByteReader._read_float(self)
        self.chirp_table_prog = ByteReader._read_int(self)
        self.lim_meas = ByteReader._read_byte(self)

        if self.lim_meas == 1:
            self.end_of_meas = ByteReader._read_int(self)
            self.scan_type = ByteReader._read_byte(self)
            self.scan_mode = ByteReader._read_byte(self)
            self.scan_dur = ByteReader._read_int(self)
            self.base_fn = ByteReader._read_string(self)

        self.arch_data = ByteReader._read_byte(self)
        self.disk_full = ByteReader._read_byte(self)
        self.sup_pow_lev = ByteReader._read_byte(self)
        self.meas_trig = ByteReader._read_byte(self)
        self.meas_start = ByteReader._read_int(self)
        self.meas_class = ByteReader._read_byte(self)
        self.store_lv0 = ByteReader._read_byte(self)
        self.store_lv1 = ByteReader._read_byte(self)

        self.rep_idx = ByteReader._read_int(self)
        self.mdf_idx = ByteReader._read_int(self)
        self.rep_num = ByteReader._read_int(self)
        self.mdf_num = ByteReader._read_int(self)
        self.mbf_name = ByteReader._read_string(self)

        self.mdf_list = []

        for i in range(self.mdf_num):
            self.mdf_list.append(ByteReader._read_string(self))

        # LV0 Data
        self.samp_t = ByteReader._read_unsigned_int(self)
        self.samp_ms = ByteReader._read_int(self)
        self.qf = ByteReader._read_byte(self)

        self.rr = ByteReader._read_float(self)
        self.rel_hum = ByteReader._read_float(self)
        self.env_t = ByteReader._read_float(self)
        self.baro_p = ByteReader._read_float(self)
        self.ws = ByteReader._read_float(self)
        self.wd = ByteReader._read_float(self)
        self.dd_volt = ByteReader._read_float(self)
        self.dd_tb = ByteReader._read_float(self)
        self.lwp = ByteReader._read_float(self)
        self.pow_if = ByteReader._read_float(self)
        self.elv = ByteReader._read_float(self)
        self.azm = ByteReader._read_float(self)
        self.status = ByteReader._read_float(self)
        self.t_pow = ByteReader._read_float(self)
        self.t_temp = ByteReader._read_float(self)
        self.r_temp = ByteReader._read_float(self)
        self.pc_temp = ByteReader._read_float(self)
        self.sky_tb = ByteReader._read_float(self)
        self.inc_el = ByteReader._read_float(self)
        self.inc_elax = ByteReader._read_float(self)

        self.t_prof = ByteReader._read_float_vector(self, self.taltn)
        self.ah_prof = ByteReader._read_float_vector(self, self.haltn)
        self.rh_prof = ByteReader._read_float_vector(self, self.haltn)

        self.s_lev = ByteReader._read_float_vector(self, self.raltn)

        prof_msk = ByteReader._read_byte_vector(self, self.raltn)

        for msk in prof_msk:

            if msk == 0:
                continue

            block_n = ByteReader._read_byte(self)

            min_block_idx = ByteReader._read_short_int_vector(self, block_n)
            max_block_idx = ByteReader._read_short_int_vector(self, block_n)

            for i in range(block_n):

                block_width = max_block_idx[i] - min_block_idx[i] + 1

                Spec = ByteReader._read_float_vector(self, block_width)

            if self.antialia == 1:
                alias_msk = ByteReader._read_byte(self)
                min_vel = ByteReader._read_float(self)

        # LV1 data

        self.ze = self.__init_nan_vector(self.raltn)
        self.mv = self.__init_nan_vector(self.raltn)
        self.sw = self.__init_nan_vector(self.raltn)
        self.sk = self.__init_nan_vector(self.raltn)
        self.kt = self.__init_nan_vector(self.raltn)
        self.ref_rat = self.__init_nan_vector(self.raltn)
        self.corr = self.__init_nan_vector(self.raltn)
        self.phi = self.__init_nan_vector(self.raltn)
        self.ze45 = self.__init_nan_vector(self.raltn)
        self.sldr = self.__init_nan_vector(self.raltn)
        self.scorr = self.__init_nan_vector(self.raltn)
        self.kdp = self.__init_nan_vector(self.raltn)
        self.diff_at = self.__init_nan_vector(self.raltn)

        for i in range(self.raltn):

            if prof_msk[i] == 0:
                continue

            self.ze[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.mv[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.sw[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.sk[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.kt[i] = float(ByteReader._read_long_long(self)) / 10**7

            if self.pol_mode == 0:
                continue

            self.ref_rat[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.corr[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.phi[i] = float(ByteReader._read_long_long(self)) / 10**7

            if self.pol_mode != 2:
                continue

            self.ze45[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.sldr[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.scorr[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.kdp[i] = float(ByteReader._read_long_long(self)) / 10**7
            self.diff_at[i] = float(ByteReader._read_long_long(self)) / 10**7

    def __init_nan_vector(self, n):

        v = np.empty(n)
        v[:] = np.nan
        return v


class Client(Status):

    def __init__(self, HostIP, Port, Password, SuppressOutput=False):

        self.__HostIP = ""
        self.__Port = 0
        self.__PWC = 0

        self.__suppressoutput = SuppressOutput

        try:
            socket.inet_aton(HostIP)
        except socket.error:
            raise Exception('Not valid IP4 address')

        if not isinstance(Port, int):
            raise Exception('The Port input variable must be an integer')

        if Port < 0 or Port > 65535:
            raise Exception(
                'The HostIP input variable must be in the range from 0 to 65535')

        self.__HostIP = HostIP
        self.__Port = Port
        self.__PWC = self.__get_password_code(Password)

    def __get_password_code(self, password):

        passcode = 0
        password = password.upper()

        for i in range(len(password)):
            passcode = passcode + ord(password[i]) * (i+1)**2

        return (passcode)

    def __send_receive(self, msgcode, filename=None, length=None, content=None):

        MESSAGE = self.__form_request_message(
            msgcode, filename=filename, length=length, content=content)

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:

            try:
                s.connect((self.__HostIP, self.__Port))
            except socket.error:
                raise Exception('Cannot connect to the server')

            s.sendall(MESSAGE)

            data = bytearray()

            while True:
                packet = s.recv(1024)
                if not packet:
                    break
                data.extend(packet)

            if self.__check_first_byte_of_response(data, msgcode) == False:
                return None

            return data

    def __check_first_byte_of_response(self, response, msgcode):

        if response[0] == 253:
            if self.__suppressoutput == False:
                print('The requested command is not known')
            return False

        if response[0] == 255:
            if self.__suppressoutput == False:
                print('The password is wrong')
            return False

        if response[0] is not msgcode:
            if self.__suppressoutput == False:
                print('Invalid response')
            return False

        return True

    def __form_request_message(self, msgcode, filename=None, length=None, content=None):

        MESSAGE = (msgcode).to_bytes(1, byteorder='little') + \
                  (self.__PWC).to_bytes(4, byteorder='little')

        if filename is not None:
            MESSAGE += bytearray(filename + chr(0), 'utf-8')

        if length is not None:
            MESSAGE += (length).to_bytes(4, byteorder='little') + \
                bytearray(content)

        return MESSAGE

    def get_radar_status(self):

        try:
            # Form the request message, 173 is the command code from the software manual
            d = self.__send_receive(173)
            S = Status(d, self.__suppressoutput)
            return S
        except:
            print('Cannot get radar status')
            return None

    def start_radar_measurements(self, mdf_name):

        try:
            # Form the request message, 171 is the command code from the software manual
            d = self.__send_receive(171, filename=mdf_name)

            if self.__suppressoutput == False:

                if d[1] == 0:
                    print('Host is not connected to radar')

                if d[1] == 1:
                    print('radar in STANDBY mode: starting measurement')

                if d[1] == 2:
                    print('radar not in STANDBY mode')

                if d[1] == 3:
                    print('Specified MDF/MBF path not valid')

                if d[1] == 4:
                    print(
                        'Host is connected to Slave radar. Slave cannot start measurement')

            return d[1]
        except:
            print('Failed to start radar measurements')
            return None

    def terminate_radar_measurements(self):

        try:
            # Form the request message, 170 is the command code from the software manual
            d = self.__send_receive(170)

            if self.__suppressoutput == False:

                if d[1] == 0:
                    print('Host is not connected to radar')

                if d[1] == 1:
                    print('no measurement running, STANDBY mode')

                if d[1] == 2:
                    print('running measurement will be terminated')

                if d[1] == 3:
                    print('cannot terminate: zero calibration running')

                if d[1] == 4:
                    print('no measurement: absolute calibration running')

                if d[1] == 5:
                    print('cannot terminate: transmitter power calibration running')

                if d[1] == 6:
                    print(
                        'Host is connected to Slave radar! Slave cannot terminate measurement')

            return d[1]
        except:
            print('Failed to terminate radar measurements')
            return None

    def start_radar_measurements_local_mdf(self, filename):

        try:
            if not isinstance(filename, str):
                print('mdf_name must be a string')

            if not Path(filename).is_file():
                print('MDF file is not found')

            f = open(filename, 'r')
            read_data = f.read()
            m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            f.close()

            d = self.__send_receive(
                172, filename=filename, length=len(read_data), content=m)

            if d[1] == 0:
                print('Host is not connected to radar')

            if d[1] == 1:
                print('radar in STANDBY mode: starting measurement')

            if d[1] == 2:
                print('radar not in STANDBY mode')

            if d[1] == 3:
                print('Host is connected to Slave radar. Slave cannot start measurement')

            return d[1]
        except:
            print('Failed to start radar measurements from local mdf')
            return None

    def get_mdf_list(self):

        try:
            d = self.__send_receive(174)
            R = MDFList(d, self.__suppressoutput)
            return R
        except:
            print('Failed to get MDF list')
            return None

    def get_radar_id(self):

        try:
            d = self.__send_receive(175)
            R = RadarID(d, self.__suppressoutput)
            return R
        except:
            print('Failed to get radar id')
            return None

    def install_local_mdf(self, filename):

        try:
            if not isinstance(filename, str):
                print('mdf_name must be a string')

            if not Path(filename).is_file():
                print('MDF file is not found')

            f = open(filename, 'r')
            read_data = f.read()
            m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            f.close()

            a = filename.split('\\')

            d = self.__send_receive(
                176, filename=a[-1] + chr(0), length=len(read_data), content=m)

            if d is None:
                return False

            if self.__suppressoutput == False:
                print('The file has been successfully transferred')

            return True
        except:
            print('Failed to install mdf on the radar')
            return False

    def get_last_sample(self):

        d = self.__send_receive(177)
        S = LastSample(d, self.__suppressoutput)
        return S


class Scan:

    def __init__(self,
                 elv=90,
                 azm=0,
                 elv_target=None,
                 azm_target=None,
                 elv_spd=1.0,
                 azm_spd=1.0,
                 align_with_wind=0,
                 ):

        if elv < 0.0 or elv > 180.0:

            raise Exception('Wrong elevation')

        if azm < 0.0 or azm > 360.0:

            raise Exception('Wrong azimuth')

        if elv_target is not None:

            if elv_target < 0.0 or elv_target > 180.0:

                raise Exception('Wrong target elevation')

        if azm_target is not None:

            if azm_target < 0.0 or azm_target > 360.0:

                raise Exception('Wrong target azimuth')

        if elv_spd < 0.0 or elv_spd > 5.0:

            raise Exception('Wrong elevation speed')

        if azm_spd < 0.0 or azm_spd > 5.0:

            raise Exception('Wrong azimuth speed')

        # constant angle
        if elv_target is None and azm_target is None:

            self.ScanType = 0
            self.ConstEl = elv
            self.ConstAz = azm

            return

        if elv_target is None or azm_target is None:

            raise Exception('Wrong targets')

        self.ScanType = 1
        self.ScanMode = 0

        self.ScanStartEl = elv
        self.ScanStopEl = elv_target
        self.ScanIncEl = 0
        self.ScanSpeedEl = elv_spd

        self.ScanStartAz = azm
        self.ScanStopAz = azm_target
        self.ScanIncAz = 0
        self.ScanSpeedAz = azm_spd

        self.AzWindAlign = align_with_wind


class MeasDefFile:

    def __init__(self):

        self.__filecode = 48856

    def __write_to_file(self, file_id, var_type, val_list):

        A = array(var_type, val_list)
        A.tofile(file_id)

    def create(self,
               filename,
               chirp_prg,
               scan_or_list_of_scans,
               # the default of None makes a frame that just goes through
               # all the defines scans once without repetion
               frames=None,
               cal_int=3600,
               timing=False,
               duration=3600,
               filelen=3600,
               # comment spirrobe: start = False (=0) means immediately
               start=False,
               start_time=0,
               # unclear how to support the ignore date, ignore hour
               trig=0,
               LV0=True,
               LV1=True,
               NoiseT=6.0,
               windfunc=4,
               basename='rpgfmcwscan_'):

        SpecCmp = True

        _windfuncs = {'Rectangle': 0,
                      'Parzen': 1,
                      'Blackman': 2,
                      'Welch': 3,
                      'Slepian2': 4,
                      'Slepian3': 5,
                      }

        if isinstance(windfunc, str):
            windfunc = _windfuncs.get(windfunc.lower(), False)
            if windfunc is False or windfunc not in list(range(6)):
                raise ValueError('Windowfunction has to be either integer or',
                                 f'a known name ({_windfuncs.keys()})')

        if isinstance(scan_or_list_of_scans, list):
            scans = scan_or_list_of_scans
        else:
            scans = [scan_or_list_of_scans]

        # Comment spirrobe:
        # assert that all scans are of the same type, albeit it is unclear
        # if this would not be supported in principle but made impossible
        # by the MDF structure
        if all(scans[0].ScanType == scan.ScanType for scan in scans):
            scantype = scans[0].ScanType
            if (scantype == 1
                    and all(scans[0].ScanMode == scan.ScanMode for scan in scans)):
                scanmode = scans[0].ScanMode
            else:
                raise ValueError("Mixing of different scanmodes",
                                 " is not possible with the MDF file structure")
        else:
            raise ValueError("Mixing of different scantypes",
                             " is not possible with the MDF file structure")

        # the framenumbers, in order
        # the framecount (default 1) as the len of the nested list
        # the start scans (default 0)
        # the stop scans (default 0)
        # the repetitions per frame (default 1)
        if frames is None:
            frames = [[0, len(scans), 1]]
        framecount = len(frames)
        framestarts = [frame[0] for frame in frames]
        # warn that it never includes the scan 0 indicating an issue
        if min(framestarts) > 0:
            print('Scan numbering starts at index 0, which was not found',
                  'in the frames definition, are you sure you defined it right?')

        framestops = [frame[1] for frame in frames]
        # warn that it goes over the number of defined scans
        if min(framestarts) > 0:
            print('Scan numbering exceeds the number of available scans',
                  ', are you sure you defined it right?')

        # ensure MDF file is always uppercase like in GUI
        filename = filename.replace(os.path.basename(filename),
                                    os.path.basename(filename).upper()
                                    )
        try:
            output_file = open(filename, 'wb')
        except OSError:
            print('Could not open file location, check if path exists')
            return

        self.__write_to_file(output_file, 'i', [
                             self.__filecode, chirp_prg, cal_int])
        self.__write_to_file(output_file, 'b', [LV0, LV1, SpecCmp, 0, 0, 0, 0])
        self.__write_to_file(output_file, 'f', [NoiseT])

        # the passing in of scans makes the most sense as a list
        # but for the MDF format the parameters of the scans are in parameter
        # order so we have to adjust this accordingly here.
        self.__write_to_file(output_file,
                             'b',
                             [scantype])
        # scan.ScanType for scan in scans])
        if scantype == 0:
            self.__write_to_file(output_file, 'f',
                                 [scan.ConstEl for scan in scans])

            self.__write_to_file(output_file, 'f',
                                 [scan.ConstAz for scan in scans])
        elif scantype == 1:
            # header type for scantype 1, scanmode refers to continuous
            self.__write_to_file(output_file, 'b', [scanmode])
            # this is the ScanCnt
            self.__write_to_file(output_file, 'i', [len(scans)])

            # elevation related 4 params
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanStartEl for scan in scans])
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanStopEl for scan in scans])
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanIncEl for scan in scans])
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanSpeedEl for scan in scans])

            # azimuth related 4 params
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanStartAz for scan in scans])
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanStopAz for scan in scans])
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanIncAz for scan in scans])
            self.__write_to_file(output_file, 'f',
                                 [scan.ScanSpeedAz for scan in scans])

            # elevation related 4 params
            self.__write_to_file(output_file, 'b',
                                 [scan.AzWindAlign for scan in scans])

            self.__write_to_file(output_file, 'i', [framecount])
            self.__write_to_file(output_file, 'i', framestarts)
            self.__write_to_file(output_file, 'i', framestops)
            self.__write_to_file(output_file, 'i', [
                                 frame[2] for frame in frames])

        # general scans must be limited in time, this is
        if scantype == 1 and timing is True:
            timing = False

        self.__write_to_file(output_file, 'b', [timing])

        # timing can either be False or True
        if timing == 0:
            self.__write_to_file(output_file, 'i', [duration])
            # add one byte to account for the null terminated byte
            self.__write_to_file(output_file, 'i', [len(basename)+1])
            # add the null terminated byte
            self.__write_to_file(output_file, 'b',
                                 [ord(i) for i in basename] + [0],
                                 )
        elif timing == 1:
            self.__write_to_file(output_file, 'i', [filelen])

        self.__write_to_file(output_file, 'b', [start])

        if start == 1:
            self.__write_to_file(output_file, 'i', [start_time, trig])

        self.__write_to_file(output_file, 'i', [windfunc, 1000])

        output_file.close()

    def output(self):

        if self.filecode != self.__filecode:
            print("No loaded MDF file")
            return
        print('#'*20, '\nGeneral parameters', '\n'+'#'*20)
        print("Chirp generator Program Index: ", self.CGProgNo)
        print("Calibration interval [s]: ", self.ZeroCallInt)
        print("Store LV0: ", self.Lev0Ena)
        print("Store LV1: ", self.Lev1Ena)
        print("Spectral compression: ", self.SpecCompEna)
        print("Store polarimetric spectral parameters in LV0: ", self.SpecParLV0Ena)
        print("File Backup: ", self.FileBackupEna)
        print("Antialiasing: ", self.AntiAliasEna)
        print("Suppress Power Leveling: ", self.PowLevSupEna)
        print("Noise Threshold Factor: ", self.NoiseThresh)

        print('\n'+'#'*20, '\nOperational parameters', '\n'+'#'*20)
        if self.ScanType == 0:
            print('Radar operates at elevation {0:.1f} and azimuth {0:.1f} deg'.format(
                self.ConstEl, self.ConstAz))

        elif self.ScanType == 1:
            print("Radar performs a general scan")

            if self.ScanMode == 0:
                print("Scan Mode: Continuous ")

            if self.ScanMode == 1:
                print("Scan Mode: Discrete, stop at each sample ")

            for i in range(self.ScanCnt):

                print('-'*20, '\nScan index: ', i)
                print('Elevation from {0:.1f} to {1:.1f} deg'.format(
                    self.ScanStartEl[i], self.ScanStopEl[i]))
                print('Elevation increment angle [deg]: ', self.ScanIncEl[i])
                print('Elevation speed [deg/s]: ', self.ScanSpeedEl[i])
                print('Azimuth from {0:.1f} to {1:.1f} deg'.format(
                    self.ScanStartAz[i], self.ScanStopAz[i]))
                print('Azimuth increment angle [deg]: ', self.ScanIncAz[i])
                print('Azimuth speed [deg/s]: ', self.ScanSpeedAz[i])
                print('Align to wind [deg/s]: ', self.AzWindAlign[i])

            print('\n'+'#'*20, '\nRepetion and order of scans', '\n'+'#'*20)
            for i in range(self.FrameCnt):

                print('Frame index: ', i)
                print('Frame start scan: ', self.FrameStartScn[i])
                print('Frame stop scan: ', self.FrameStopScn[i])
                print('Frame repetition number: ', self.FrameRep[i])

        if self.Timing == 0:

            print('Limited measurement')
            print('Duration [s]: ', self.Duration)

        if self.Timing == 1:

            print('Unlimited measurement')
            print('File Length [s]', self.FileLen)

        if self.MeasStart == 1:

            print('Start time: ', self.StartTime)
            print('Trigger condition: ', self.TrigCond)

        if self.WindFunc == 0:

            print('Window Function: Rectangle')

        if self.WindFunc == 1:

            print('Window Function: Parzen')

        if self.WindFunc == 2:

            print('Window Function: Blackman')

        if self.WindFunc == 3:

            print('Window Function: Welch')

        if self.WindFunc == 4:

            print('Window Function: Slepian2')

        if self.WindFunc == 5:

            print('Window Function: Slepian3')

        print('ADC voltage range [mV]: ', self.ADCVoltRng/1000)

    def read(self, filepath):

        with open(filepath, mode='rb') as file:
            # can close the file here already ..
            data = file.read()

        R = ByteReader(data, 0)

        self.filecode = R._read_int()

        if self.filecode != 48856:
            raise Exception("Filecode does not correspond to MDF format")

        self.CGProgNo = R._read_int()
        self.ZeroCallInt = R._read_int()
        self.Lev0Ena = R._read_byte()
        self.Lev1Ena = R._read_byte()
        self.SpecCompEna = R._read_byte()
        self.SpecParLV0Ena = R._read_byte()
        self.FileBackupEna = R._read_byte()
        self.AntiAliasEna = R._read_byte()
        self.PowLevSupEna = R._read_byte()
        self.NoiseThresh = R._read_float()
        self.ScanType = R._read_byte()

        if self.ScanType == 0:
            self.ConstEl = R._read_float()
            self.ConstAz = R._read_float()

        if self.ScanType == 1:
            self.ScanMode = R._read_byte()
            self.ScanCnt = R._read_int()

            self.ScanStartEl = []
            self.ScanStopEl = []
            self.ScanIncEl = []
            self.ScanSpeedEl = []
            self.ScanStartAz = []
            self.ScanStopAz = []
            self.ScanIncAz = []
            self.ScanSpeedAz = []
            self.AzWindAlign = []
            # 8 vars for each scan to be read, 4 elv, 4 az
            scandata = [[R._read_float() for j in range(self.ScanCnt)]
                        for i in range(8)]

            # elevation related params
            self.ScanStartEl.extend(scandata[0])
            self.ScanStopEl.extend(scandata[1])
            self.ScanIncEl.extend(scandata[2])
            self.ScanSpeedEl.extend(scandata[3])

            # azimuth related params
            self.ScanStartAz.extend(scandata[4])
            self.ScanStopAz.extend(scandata[5])
            self.ScanIncAz.extend(scandata[6])
            self.ScanSpeedAz.extend(scandata[7])

            # whether to align the azimuth into the wind,
            # for each scan a byte
            self.AzWindAlign.extend([R._read_byte()
                                    for j in range(self.ScanCnt)])

            # comment spirrobe
            # maybe this needs to be refactored to in case there are
            # several frames. requires testing.
            self.FrameCnt = R._read_int()
            self.FrameStartScn = []
            self.FrameStopScn = []
            self.FrameRep = []

            self.FrameStartScn.extend([R._read_int()
                                       for j in range(self.FrameCnt)])
            self.FrameStopScn.extend([R._read_int()
                                      for j in range(self.FrameCnt)])
            self.FrameRep.extend([R._read_int()
                                  for j in range(self.FrameCnt)])

        self.Timing = R._read_byte()
        if self.Timing == 0:
            self.Duration = R._read_int()
            self.BaseNmLen = R._read_int()
            self.BaseNm = ''.join([chr(R._read_byte())
                                   for j in range(self.BaseNmLen)])

        elif self.Timing == 1:
            self.FileLen = R._read_int()

        self.MeasStart = R._read_byte()
        if self.MeasStart == 1:

            self.StartTime = R._read_int()
            self.TrigCond = R._read_int()

        self.WindFunc = R._read_int()
        # comment spirrobe:
        # is this really the ADCVoltRng in an RPG FMCW?
        # here in the RPG FMCW it is always 1000 (see also .create method)
        self.ADCVoltRng = R._read_int()

        if R._ByteReader__i == R._ByteReader__input.__len__():
            print(f'Finished processing {filepath}')
        else:
            # this is pure guesswork based on making different MDF files via
            # the gui
            self.extra = R._read_signed_short()

            # comment spirrobe:
            # this was the number I found during testing, as I do not know
            # what the 256 stands for in this case there could be other cases
            # where something else is defined in the MDF.
            if self.extra == 256:
                self.BaseNmLen, self.BaseNm = R._read_int(), R._read_string()

            print(f'{R._ByteReader__i} bytes of',
                  f'{R._ByteReader__input.__len__()} bytes read')
            print(f'{R._ByteReader__input[R._ByteReader__i:]} remains')


class MeasBatchFile:

    def __init__(self):

        self.__filecode = 58856

    def __write_to_file(self, file_id, var_type, val_list):

        A = array(var_type, val_list)
        A.tofile(file_id)

    def create(self,
               filename,
               mdflist,
               repetitions=1,
               ):

        if isinstance(mdflist, list):
            mdflist = mdflist
        else:
            mdflist = [mdflist]

        # ensure MBF filename is always uppercase
        filename = filename.replace(os.path.basename(filename),
                                    os.path.basename(filename).upper()
                                    )
        try:
            output_file = open(filename, 'wb')
        except OSError:
            print('Could not open file location, check if path exists')
            return

        self.__write_to_file(output_file, 'i', [self.__filecode,])
        self.__write_to_file(output_file, 'i', [len(mdflist)])

        for mdffile in mdflist:
            shortname = mdffile.replace('\\', '/').split('/')[-1]
            fullname = mdffile.replace('/', '\\')
            # add one byte to account for the null terminated byte
            self.__write_to_file(output_file, 'i', [len(fullname)+1])
            # add the null terminated byte
            self.__write_to_file(output_file, 'b',
                                 [ord(i) for i in fullname] + [0],
                                 )
            # add one byte to account for the null terminated byte
            self.__write_to_file(output_file, 'i', [len(shortname)+1])
            # add the null terminated byte
            print(shortname, [len(shortname)+1])
            self.__write_to_file(output_file, 'b',
                                 [ord(i) for i in shortname] + [0],
                                 )
        # self.__write_to_file(output_file, 'b',
        #                      [0],
        #                      )
        self.__write_to_file(output_file, 'i', [repetitions,])
        output_file.close()

    def output(self):

        if self.filecode != self.__filecode:
            print("No loaded MBF file")
            return

        print('#'*20, '\nList of MDFs', '\n'+'#'*20)

        for i in range(self.nMDFs):
            print(f'MDF file {i+1}', self.MDFlist[i])

        print(f'Batch will be repeated: {self.Repetitions}')

    def read(self, filepath):

        with open(filepath, mode='rb') as file:
            # can close the file here already ..
            data = file.read()

        R = ByteReader(data, 0)

        self.filecode = R._read_int()

        if self.filecode != 58856:
            raise Exception("Filecode does not correspond to MBF format")

        self.nMDFs = R._read_int()
        print(self.nMDFs)
        if self.nMDFs:
            self.MDFlist = []

        for i in range(self.nMDFs):
            FullNameLen = R._read_int()
            FullName = ''.join([chr(R._read_byte())
                                for j in range(FullNameLen)][:-1])

            ShortNameLen = R._read_int()
            ShortName = ''.join([chr(R._read_byte())
                                 for j in range(ShortNameLen)][:-1])
            self.MDFlist += [[FullName, ShortName]]

        self.Repetitions = R._read_int()

        if R._ByteReader__i == R._ByteReader__input.__len__():
            print(f'Finished processing {filepath}')
        else:
            print(f'{R._ByteReader__i} bytes of',
                  f'{R._ByteReader__input.__len__()} bytes read')
            print(f'{R._ByteReader__input[R._ByteReader__i:]} remains')


def get_radar_status(HOSTIP, PORT, PASSWORD):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.get_radar_status()


def start_radar_measurements(HOSTIP, PORT, PASSWORD, MDF_FILENAME):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.start_radar_measurements(MDF_FILENAME)


def start_radar_measurements_local_mdf(HOSTIP, PORT, PASSWORD, MDF_FILENAME):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.start_radar_measurements_local_mdf(MDF_FILENAME)


def terminate_radar_measurements(HOSTIP, PORT, PASSWORD):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.terminate_radar_measurements()


def get_mdf_list(HOSTIP, PORT, PASSWORD):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.get_mdf_list()


def get_radar_id(HOSTIP, PORT, PASSWORD):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.get_radar_id()


def install_local_mdf(HOSTIP, PORT, PASSWORD, MDF_FILENAME):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD)
    X.install_local_mdf(MDF_FILENAME)


def get_last_sample(HOSTIP, PORT, PASSWORD):

    PORT = int(PORT)

    X = Client(HOSTIP, PORT, PASSWORD, SuppressOutput=True)
    X.get_last_sample()


if __name__ == '__main__':

    status = False

    if len(sys.argv) > 1:
        if sys.argv[1] == "get_radar_status":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print(
                        'python RadarControl.py get_radar_status HOSTIP PORT PASSWORD')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    status = True

        if sys.argv[1] == "start_radar_measurements":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print(
                        'python RadarControl.py start_radar_measurements HOSTIP PORT PASSWORD MDF_FILENAME')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    print('MDF_FILENAME is a string with a MDF/MBF filename on the host PC in the format *.MDF/*.MBF. The string can also contain the full path if the needed file is not in the default MDF / MBF directory.')
                    print(
                        'In order to get the list of available MDF/MBF files in the default MDF / MBF directory, refer to the get_mdf_list command.')
                    status = True

        if sys.argv[1] == "start_radar_measurements_local_mdf":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print(
                        'python RadarControl.py start_radar_measurements_local_mdf HOSTIP PORT PASSWORD MDF_FILENAME')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    print('MDF_FILENAME is a string with a MDF/MBF filename on the local PC in the format path/*.MDF (path/*.MBF). The string must contain the full path. The MDF/MBF file will be transferred to the host PC and launched. Please note, that the file will NOT be stored in the default MDF/MBF folder on the host PC.')
                    status = True

        if sys.argv[1] == "terminate_radar_measurements":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print(
                        'python RadarControl.py terminate_radar_measurements HOSTIP PORT PASSWORD')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    status = True

        if sys.argv[1] == "get_mdf_list":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print('python RadarControl.py get_mdf_list HOSTIP PORT PASSWORD')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    status = True

        if sys.argv[1] == "get_radar_id":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print('python RadarControl.py get_radar_id HOSTIP PORT PASSWORD')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    status = True

        if sys.argv[1] == "install_local_mdf":
            if len(sys.argv) == 3:
                if sys.argv[2] == "help":
                    print('\nCommand template:\n')
                    print(
                        'python RadarControl.py install_local_mdf HOSTIP PORT PASSWORD MDF_FILENAME')
                    print(
                        '\nHOSTIP is the IP address of the host PC in format a.b.c.d, where a,b,c,d are numbers from 0 to 255.')
                    print(
                        'PORT is the port of the Data Server on the host PC. Must be a positive integer number.')
                    print(
                        'PASSWORD is a string with the User Password. If the user password is not set, type any symbol.')
                    print('MDF_FILENAME is a string with a MDF/MBF filename on the local PC in the format path/*.MDF (path/*.MBF). The string must contain the full path. The MDF/MBF file will be transferred and stores in the default MDF/MBF folder on the host PC. Please note, that the measurements will NOT be started. Please refer to the start_radar_measurements command.')
                    status = True

    if len(sys.argv) == 5:
        globals()[sys.argv[1]](sys.argv[2], sys.argv[3], sys.argv[4])
        status = True

    if len(sys.argv) == 6:
        globals()[sys.argv[1]](sys.argv[2],
                               sys.argv[3], sys.argv[4], sys.argv[5])
        status = True

    if status == False:
        print('\nError: Wrong input variables.\n')
        print('The command line interface of the RadarControl module supports the following commands:\n')
        print('get_radar_status')
        print('start_radar_measurements')
        print('start_radar_measurements_local_mdf')
        print('terminate_radar_measurements')
        print('get_mdf_list')
        print('get_radar_id')
        print('install_local_mdf')

        print('\nIn order to get help for a certain command use the following command:')
        print('python RadarControl.py CommandName "help"')
