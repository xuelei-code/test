# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 10:52:38 2018

@author: lei
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 09:23:18 2018

@author: lei
"""

#import matplotlib as mpl  
#mpl.use('Agg')  
import numpy as np  
import matplotlib.pyplot as plt  
import os, shutil
import csv
import math
import time
from obspy.core import UTCDateTime
from os import listdir
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.clients.filesystem import sds
from obspy import read
import pandas as pd
import concurrent.futures

# ratio、filepath和outputpath是根据实际情况需要修改的
#设定烈度计换算系数，jopens系统为50000，earthworm系统为1671.8
ratio = 1671.8
# SDS格式的数据存储位置
filepath='//10.35.176.56/public'
# 普通文件夹格式存储位置
#filepath = 'F:/work/mycode/mseedinput/2008117'
# outputpath是输出结果的存储位置
outputpath = 'F:/work/mycode/intensity_result'

#以下内容请不要随便修改
#获取开始时间
time_begin = time.time()

#把发震信息读到列表message中
message = []
for line in open("auto_eventdata"): 
    message.append(line.replace('\n',''))
    
#震源信息读出
eq_year = int(message[0])
eq_month = int(message[1])
eq_date = int(message[2])
eq_hour = int(message[3])
eq_min = int(message[4])
eq_sec = int(message[5])
eq_lat = float(message[6])
eq_lon = float(message[7])
eq_mag = message[8]
eq_dep = float(message[9])

# 发震时间为t0
t0 = UTCDateTime(str(eq_year) + '-'+  str(eq_month) + '-'+  str(eq_date) + 'T' + str(eq_hour) + ':' + str(eq_min) + ':' + str(eq_sec))
eq_jday = t0.julday

# 根据station.csv中的信息把台站名,经度，纬度，台站地址取出
with open("stationdata.csv", "r", encoding = "utf-8") as f1:
    reader = csv.reader(f1)
    sta_name = [row[0] for row in reader]
with open("stationdata.csv", "r", encoding = "utf-8") as f2:
    reader = csv.reader(f2)
    sta_lat = [row[1] for row in reader]
with open("stationdata.csv", "r", encoding = "utf-8") as f3:
    reader = csv.reader(f3)
    sta_lon = [row[2] for row in reader]

# 读入波形数据
#SDS文件系统的数据读入
#obspy中SDS读取接口读取
#path = sds.Client(filepath)
#stEs = path.get_waveforms("FJ", "*", "*", "EIE", t0-10, t0+290).sort().merge(fill_value=0, interpolation_samples=1)
#stNs = path.get_waveforms("FJ", "*", "*", "EIN", t0-10, t0+290).sort().merge(fill_value=0, interpolation_samples=1)
#stZs = path.get_waveforms("FJ", "*", "*", "EIZ", t0-10, t0+290).sort().merge(fill_value=0, interpolation_samples=1)


# 正常文件夹中的mseed数据读入
#sts = read(filepath + "/*.mseed",starttime=t0-10,endtime=t0+290).sort()

#建立文件索引读取
mseedname_E = []
mseedname_N = []
mseedname_Z = []
for name in sta_name:
    sta_name_E = filepath + "/" + str(eq_year) + "/FJ/" + name + "/EIE.D/" + "FJ." + name + ".00.EIE.D." + str(eq_year) + "." + str(eq_jday)
    sta_name_N = filepath + "/" + str(eq_year) + "/FJ/" + name + "/EIN.D/" + "FJ." + name + ".00.EIN.D." + str(eq_year) + "." + str(eq_jday)
    sta_name_Z = filepath + "/" + str(eq_year) + "/FJ/" + name + "/EIZ.D/" + "FJ." + name + ".00.EIZ.D." + str(eq_year) + "." + str(eq_jday)
    if os.path.exists(sta_name_E) and os.path.exists(sta_name_N) and os.path.exists(sta_name_Z):
        mseedname_E.append(sta_name_E)
        mseedname_N.append(sta_name_N)
        mseedname_Z.append(sta_name_Z)
        
def read_sdsfile(mseedfilename):
    sts = read(mseedfilename, starttime=t0-10,endtime=t0+290)
    return sts
sts_E = read()
sts_N = read()
sts_Z = read()
for name_E in mseedname_E:
    sts_E += read_sdsfile(name_E)
for name_N in mseedname_N:
    sts_N += read_sdsfile(name_N)
for name_Z in mseedname_Z:
    sts_Z += read_sdsfile(name_Z)
    
sts_E.sort().merge(fill_value=0, interpolation_samples=1)
sts_N.sort().merge(fill_value=0, interpolation_samples=1)
sts_Z.sort().merge(fill_value=0, interpolation_samples=1)

time2 = time.time()
print("数据读取完毕！花费时间为%d秒!"  % (time2 - time_begin))