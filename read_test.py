# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:46:19 2018

@author: lei
"""

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

# ratio、filepath和outputpath是根据实际情况需要修改的
#设定烈度计换算系数，jopens系统为50000，earthworm系统为1671.8
ratio = 1671.8
# SDS格式的数据存储位置
#filepath='//10.35.176.56/public'
# 普通文件夹格式存储位置
filepath = '//10.35.176.56/public'
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
print(t0)

# mseed文件读入
# 遇到发震时间为跨天的情况
# 当发震时间为这一天的00：00-00：10分之间时
t1 = t0 - 10
t1_year = t1.year
t1_jday = t1.julday
t2 = t0 + 290
t2_year = t2.year
t2_jday = t2.julday
#if eq_hour == 0 and eq_min == 0 and eq_sec <= 10:
#    sts = read(filepath + "/" + str(t1_year) + "/" + str(t1_jday) + "/*",starttime=t1)
#    sts += read(filepath + "/" + str(eq_year) + "/" + str(eq_jday) + "/*",endtime=t0 + 290)
## 当发震时间为这一天23：55：10以后时，t2为读取结束的时间点    
#elif (eq_hour == 23 and eq_min > 55) or (eq_hour == 23 and eq_min == 55 and eq_sec >10):
#    sts = read(filepath + "/" + str(eq_year) + "/" + str(eq_jday) + "/*",starttime=t1)
#    sts += read(filepath + "/" + str(t2_year) + "/" + str(t2_jday) + "/*",endtime=t2)
## 没有以上特殊情况时
#else:
#    sts = read(filepath + "/" + str(eq_year) + "/" + str(eq_jday) + "/*",starttime=t0-10,endtime=t0+290)

# 正常情况下
if t1_jday == t2_jday:
    sts = read(filepath + "/" + str(eq_year) + "/" + str(eq_jday) + "/*",starttime=t1,endtime=t2)
# 如果有跨天的情况
else:
    sts = read(filepath + "/" + str(t1_year) + "/" + str(t1_jday) + "/*",starttime=t1)
    sts += read(filepath + "/" + str(t2_year) + "/" + str(t2_jday) + "/*",endtime=t2)

stZs = sts.select(channel = "EIZ").merge(fill_value=0, interpolation_samples=1).sort()

time1 = time.time()
print("数据读取完毕！花费时间为%d秒!"  % (time1 - time_begin))

print(stZs)