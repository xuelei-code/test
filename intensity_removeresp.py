# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 15:05:14 2018

@author: lei
"""

import os
import re
import csv
import math
from obspy import read, read_inventory
from obspy.core import UTCDateTime
from os import listdir
from shutil import copyfile
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
import numpy as np

#把发震信息读到列表message中
message = []
for line in open("auto_eventdata"): 
    message.append(line.replace('\n',''))
    
# 发震时间为t0
# t0 = UTCDateTime(message[0] + '-'+  message[1] + '-'+  message[2] + 'T' + message[3] + ':' + message[4] + ':' + message[5])
t0 = UTCDateTime("2018-04-22T16:30:00.000")
#震源信息读出
eq_lat = float(message[6])
eq_lon = float(message[7])
eq_mag = message[8]
eq_dep = float(message[9])

# 波形数据存储位置
filepath='F:/work/mycode/test/mseed/'

# 读入仪器响应文件
inv = read_inventory('FJ.dataless.xml')
pre_filt = [0.05, 0.1, 10, 20]

# 读入波形数据,去除仪器响应,输出为速度文件
sts=read(filepath+'*',starttime = t0 - 60, endtime = t0 + 300).merge().sort()
V_sts = sts.remove_response(inventory=inv, output='VEL', pre_filt=pre_filt)
V_sts.write("VEL.mseed", format="MSEED")

# 读入波形数据,去除仪器响应,输出为加速度文件
sts=read(filepath+'*',starttime = t0 - 60, endtime = t0 + 300).merge().sort()
print(np.max(sts[0].data))
A_sts = sts.remove_response(inventory=inv, output='ACC', pre_filt=pre_filt)
A_sts.write("ACC.mseed", format="MSEED")
print(np.max(A_sts[0].data))

# 把台站名,经度，纬度取出
station_list = []
with open("stationdata.csv", "r", encoding = "utf-8") as f:
    reader = csv.reader(f)
    sta_name = [row[0] for row in reader]
with open("stationdata.csv", "r", encoding = "utf-8") as f:
    reader = csv.reader(f)
    sta_lat = [row[1] for row in reader]
with open("stationdata.csv", "r", encoding = "utf-8") as f:
    reader = csv.reader(f)
    sta_lon = [row[2] for row in reader]
    
#以下开始读入速度的mseed文件，进行pgv的计算
sts=read('VEL.mseed')
stEs = sts.select(channel = "EIE")
stNs = sts.select(channel = "EIN")
stZs = sts.select(channel = "EIZ")

# mseed文件中有的台站名为
stas=[]
for sts in stEs:
    stas.append(sts.stats.station)
    
# 给台站经纬度赋值    
for sta in stas:
    for i in range(len(sta_name)):
        if sta == sta_name[i]:
            stZs.select(station=sta)[0].stats.stla = float(sta_lat[i])
            stZs.select(station=sta)[0].stats.stlo = float(sta_lon[i])
            stEs.select(station=sta)[0].stats.stla = float(sta_lat[i])
            stEs.select(station=sta)[0].stats.stlo = float(sta_lon[i])
            stNs.select(station=sta)[0].stats.stla = float(sta_lat[i])
            stNs.select(station=sta)[0].stats.stlo = float(sta_lon[i])

# 计算震中距，并赋值    
for sts in stEs:
    sts.stats.distance = gps2dist_azimuth(sts.stats.stla,sts.stats.stlo,eq_lat,eq_lon)
for sts in stNs:
    sts.stats.distance = gps2dist_azimuth(sts.stats.stla,sts.stats.stlo,eq_lat,eq_lon)
for sts in stZs:
    sts.stats.distance = gps2dist_azimuth(sts.stats.stla,sts.stats.stlo,eq_lat,eq_lon)

# 计算理论到时，把理论到时之前10s的平均值作为矫正值
# iasp91模型 
model = TauPyModel(model="iasp91")
for sts in stEs:
    arrivals = model.get_travel_times_geo(source_depth_in_km = eq_dep, source_latitude_in_deg = eq_lat, source_longitude_in_deg = eq_lon, receiver_latitude_in_deg = sts.stats.stla, receiver_longitude_in_deg = sts.stats.stlo)
    arr = arrivals[0]
    cor_st = int((arr.time + 50)*sts.stats.sampling_rate)
    cor_en = int((arr.time + 60)*sts.stats.sampling_rate)
    cor = np.mean(sts.data[cor_st:cor_en])
    sts.data = np.array([x - cor for x in sts.data])
for sts in stNs:
    arrivals = model.get_travel_times_geo(source_depth_in_km = eq_dep, source_latitude_in_deg = eq_lat, source_longitude_in_deg = eq_lon, receiver_latitude_in_deg = sts.stats.stla, receiver_longitude_in_deg = sts.stats.stlo)
    arr = arrivals[0]
    cor_st = int((arr.time + 50)*sts.stats.sampling_rate)
    cor_en = int((arr.time + 60)*sts.stats.sampling_rate)
    cor = np.mean(sts.data[cor_st:cor_en])
    sts.data = np.array([x - cor for x in sts.data])
for sts in stZs:
    arrivals = model.get_travel_times_geo(source_depth_in_km = eq_dep, source_latitude_in_deg = eq_lat, source_longitude_in_deg = eq_lon, receiver_latitude_in_deg = sts.stats.stla, receiver_longitude_in_deg = sts.stats.stlo)
    arr = arrivals[0]
    cor_st = int((arr.time + 50)*sts.stats.sampling_rate)
    cor_en = int((arr.time + 60)*sts.stats.sampling_rate)
    cor = np.mean(sts.data[cor_st:cor_en])
    sts.data = np.array([x - cor for x in sts.data])
    
#按台站名分别计算各个方向的pgv和合成pgv
pgv_EW = []
pgv_NS = []
pgv_UD = []
pgv_thr = []
I_pgv = []
for sta in stas:
    sts_E = stEs.select(station = sta)
    sts_N = stNs.select(station = sta)
    sts_Z = stZs.select(station = sta)
#    print(sts_E[0].stats)
#    print(np.max(sts_E[0].data))
#    sts_E.plot()
    EW = np.max(np.abs(sts_E[0].data))
    NS = np.max(np.abs(sts_N[0].data))
    UD = np.max(np.abs(sts_Z[0].data))
    thr = np.max(np.sqrt(np.square(sts_E[0].data) + np.square(sts_N[0].data) + np.square(sts_Z[0].data)))
    i_pgv = 3.00*math.log10(thr) + 9.77
    pgv_EW.append(EW)
    pgv_NS.append(NS)
    pgv_UD.append(UD)
    pgv_thr.append(thr)
    I_pgv.append(i_pgv)
#print(pgv_EW)
#print(pgv_NS)
#print(pgv_UD)
#print(pgv_thr)
#print(I_pgv)