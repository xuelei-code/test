# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:16:38 2018

@author: lei
"""

import os
import re
import csv
from obspy import read, read_inventory
from obspy.core import UTCDateTime
from os import listdir
from shutil import copyfile
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
import numpy as np

filepath='F:/work/mycode/test/mseed/'
#把发震信息读到列表message中
message = []
for line in open("auto_eventdata"): 
    message.append(line.replace('\n',''))
    
eq_lat = float(message[6])
eq_lon = float(message[7])
eq_mag = message[8]
eq_dep = float(message[9])
# 读入时间
t0 = UTCDateTime("2018-04-22T16:30:00.000")
#设置文件读入路径
filepath='F:/work/mycode/test/mseed/'
# 读入仪器响应文件
inv = read_inventory('FJ.dataless.xml')
pre_filt = (0.005, 0.1, 10, 11)
# 读入波形数据,去除仪器响应,输出为速度波形
sts=read(filepath+'*',starttime = t0 - 60, endtime = t0 + 300).merge().sort()
V_sts = sts.remove_response(inventory=inv, output='VEL', pre_filt=pre_filt)
V_sts.write("VEL.mseed", format="MSEED")
# 读入波形数据,去除仪器响应,输出为加速度波形
sts=read(filepath+'*',starttime = t0 - 60, endtime = t0 + 300).merge().sort()
A_sts = sts.remove_response(inventory=inv, output='ACC', pre_filt=pre_filt)
A_sts.write("ACC.mseed", format="MSEED")

sts=read('VEL.mseed')
E_sts = sts.select(channel = "EIE")
print(E_sts[1].stats)
