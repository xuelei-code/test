# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 09:23:18 2018

@author: lei
"""

import os
import csv
import math
import numpy
import time
from obspy import read, read_inventory
from obspy.core import UTCDateTime
from os import listdir
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.clients.filesystem import sds
import numpy as np
import pandas as pd

#设定烈度计换算系数，jopens系统为50000，earthworm系统为1671.8
ratio = 1671.8
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
# t0 = UTCDateTime("2018-04-22T16:30:00.000")

# 波形数据存储位置
filepath='//10.35.176.56/public'

## 读入仪器响应文件
#inv = read_inventory('FJ.dataless.xml')
#pre_filt = [0.05, 0.1, 10, 20]

# 读入波形数据,输出为加速度文件
path = sds.Client(filepath)
sts = path.get_waveforms("FJ", "*", "*", "*", t0-60, t0+300).sort()
A_sts = sts.filter('bandpass', freqmin=1.0, freqmax=10, zerophase=True)
for i in range(len(listdir(filepath))):
    A_sts[i].data = np.array([x / ratio for x in sts[i].data])
#    print(np.max(A_sts[i].data))
A_sts.write("ACC.mseed", format="MSEED")
#将加速度文件积分，得到速度文件保存
V_sts = A_sts.integrate()
V_sts.write("VEL.mseed", format="MSEED")

# 把台站名,经度，纬度，台站地址取出
station_list = []
with open("stationdata.csv", "r", encoding = "utf-8") as f1:
    reader = csv.reader(f1)
    sta_name = [row[0] for row in reader]
with open("stationdata.csv", "r", encoding = "utf-8") as f2:
    reader = csv.reader(f2)
    sta_lat = [row[1] for row in reader]
with open("stationdata.csv", "r", encoding = "utf-8") as f3:
    reader = csv.reader(f3)
    sta_lon = [row[2] for row in reader]
    
#以下开始读入速度的mseed文件，进行pgv的计算
sts=read('VEL.mseed')
stEs = sts.select(channel = "EIE").merge().sort()
stNs = sts.select(channel = "EIN").merge().sort()
stZs = sts.select(channel = "EIZ").merge().sort()

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
    EW = np.max(np.abs(sts_E[0].data))
    NS = np.max(np.abs(sts_N[0].data))
    UD = np.max(np.abs(sts_Z[0].data))
    thr = np.max(np.sqrt(np.square(sts_E[0].data) + np.square(sts_N[0].data) + np.square(sts_Z[0].data)))
    if thr ==0:
        i_pgv = 0
        print("%s的I-pgv的值为空" % stEs[0].stats.station)
    else:
        i_pgv = 3.00*math.log10(thr) + 9.77
    pgv_EW.append(EW)
    pgv_NS.append(NS)
    pgv_UD.append(UD)
    pgv_thr.append(thr)
    I_pgv.append(i_pgv)
    
#以下开始读入加速度的mseed文件，进行pga的计算
sts=read('ACC.mseed')
stEs = sts.select(channel = "EIE").merge().sort()
stNs = sts.select(channel = "EIN").merge().sort()
stZs = sts.select(channel = "EIZ").merge().sort()

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
    
#按台站名分别计算各个方向的pga和合成pga
pga_EW = []
pga_NS = []
pga_UD = []
pga_thr = []
I_pga = []
for sta in stas:
    sts_E = stEs.select(station = sta)
    sts_N = stNs.select(station = sta)
    sts_Z = stZs.select(station = sta)
    EW = np.max(np.abs(sts_E[0].data))
    NS = np.max(np.abs(sts_N[0].data))
    UD = np.max(np.abs(sts_Z[0].data))
    thr = np.max(np.sqrt(np.square(sts_E[0].data) + np.square(sts_N[0].data) + np.square(sts_Z[0].data)))
    if thr == 0:
        i_pga = 0
        print("%s的I-pga的值为空" % stEs[0].stats.station)
    else:
        i_pga = 3.17*math.log10(thr) + 6.59
    pga_EW.append(EW)
    pga_NS.append(NS)
    pga_UD.append(UD)
    pga_thr.append(thr)
    I_pga.append(i_pga)
    
#取出经纬度，震中距,台站名等信息
st_lat = []
st_lon = []
distance = []
for sta in stas:
    st_lat.append(stEs.select(station = sta)[0].stats.stla)
    st_lon.append(stEs.select(station = sta)[0].stats.stlo)
    distance.append((stEs.select(station = sta)[0].stats.distance[0])/1000)
    
#根据I_pgv和I_pga的值对烈度I进行修正
I = []
for i in range(len(stas)):
    if I_pga[i] >= 6.0 and I_pgv[i] >=6.0:
        I.append(I_pgv[i])
    else:
        I.append((I_pga[i] + I_pgv[i]) / 2)
for i in range(len(stas)):
    if I[i] < 1.0:
        I[i] = 1.0
    if I[i] > 12.0:
        I[i] = 12.0
        
# 把结果写入到文件
r1 = np.array([stas, st_lat, st_lon, distance, pga_EW, pga_NS, pga_UD, pga_thr, pgv_EW, pgv_NS, pgv_UD, pgv_thr, I])
r1 = numpy.ndarray.tolist(r1.T)
result = sorted(r1, key=lambda r1_tuple: r1_tuple[3])
name = ['station', 'lat', 'lon', 'distance(km)', 'EW-PGA', 'NS-PGA', 'UD-PGA', 'PGA', 'EW-PGV', 'NS-PGV', 'UD-PGV', 'PGV', 'I']
intensity = pd.DataFrame(columns=name, data=result)
intensity.to_csv('py_intensity.csv')

#输出intensitydata和epicenterdata用于gmt作图
with open('py_intensity.csv','r') as r:
    lines=r.readlines()
with open('intensitydata','w') as w: 
    w.writelines(lines[1:])
    
with open('intensitydata','r') as r:
    lines=r.readlines()
with open('intensitydata','w',encoding="utf-8") as f_w:
    for line in lines:
        if "," in line:
            line = line.replace(",", " ")
        f_w.write(line)
        
with open('epicenterdata','w') as w:
    w.write(str(eq_year) + '/'+  str(eq_month) + '/'+  str(eq_date) + ' ' + 'M' + str(eq_mag) +' ' + str(eq_lat) + ' ' + str(eq_lon) + ' ' + str(eq_dep))

#调用gmt画图
os.system('./pga.gmt intensitydata epicenterdata')
os.system('./pgv.gmt intensitydata epicenterdata')
os.system('./intensity.gmt intensitydata epicenterdata')
os.remove('intensitydata')
os.remove('epicenterdata')
os.remove('ACC.mseed')
os.remove('VEL.mseed')

# 获取结束时间
time_over = time.time()
print("计算完毕！花费时间为%d秒!"  % (time_over - time_begin))
