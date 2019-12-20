#!/usr/bin/env python
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

# 创建数据输出文件夹
dirname = message[0] + message[1] + message[2] + 'M' + message[8]
dirnamelist = listdir(outputpath)
i = 0
for name in dirnamelist:
    if dirname in name or dirname == name:
        i = i + 1
if i > 0:
    dirname = dirname + '_' + str(i)
outputdir = outputpath + '/' + dirname
os.mkdir(outputdir)

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
with open("stationdata.csv", "r", encoding = "utf-8") as f4:
    reader = csv.reader(f4)
    sta_loc= [row[4] for row in reader]

# 读入波形数据
# SDS文件系统的数据读入
#path = sds.Client(filepath)
#sts = path.get_waveforms("FJ", "*", "*", "*", t0-10, t0+290).sort()
    
# mseed文件读入

t1 = t0 - 10
t1_year = t1.year
t1_jday = t1.julday
t2 = t0 + 290
t2_year = t2.year
t2_jday = t2.julday

# 正常情况下
if t1_jday == t2_jday:
    sts = read(filepath + "/" + str(eq_year) + "/" + str(eq_jday) + "/*",starttime=t1,endtime=t2)
# 如果有跨天的情况
else:
    sts = read(filepath + "/" + str(t1_year) + "/" + str(t1_jday) + "/*",starttime=t1)
    sts += read(filepath + "/" + str(t2_year) + "/" + str(t2_jday) + "/*",endtime=t2)

time1 = time.time()
print("数据读取完毕！花费时间为%d秒!"  % (time1 - time_begin))

# 以下开始进行pga的计算
# 读入各方向上的波形，滤波
stEs = sts.select(channel = "*E").merge(fill_value=0, interpolation_samples=1).sort()
stNs = sts.select(channel = "*N").merge(fill_value=0, interpolation_samples=1).sort()
stZs = sts.select(channel = "*Z").merge(fill_value=0, interpolation_samples=1).sort()

for ssts in stEs:
    if hasattr(ssts.data,'mask'):
        ssts.data=ssts.data.filled(0)
for ssts in stNs:
    if hasattr(ssts.data,'mask'):
        ssts.data=ssts.data.filled(0)
for ssts in stZs:
    if hasattr(ssts.data,'mask'):
        ssts.data=ssts.data.filled(0) 
stEs.filter('bandpass', freqmin=1.0, freqmax=10)
stNs.filter('bandpass', freqmin=1.0, freqmax=10)
stZs.filter('bandpass', freqmin=1.0, freqmax=10) 
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
            stZs.select(station=sta)[0].stats.loc = sta_loc[i]
            stEs.select(station=sta)[0].stats.stla = float(sta_lat[i])
            stEs.select(station=sta)[0].stats.stlo = float(sta_lon[i])
            stEs.select(station=sta)[0].stats.loc = sta_loc[i]
            stNs.select(station=sta)[0].stats.stla = float(sta_lat[i])
            stNs.select(station=sta)[0].stats.stlo = float(sta_lon[i])
            stNs.select(station=sta)[0].stats.loc = sta_loc[i]

# 计算震中距，并赋值    
for sts in stEs:
    sts.stats.distance = gps2dist_azimuth(sts.stats.stla,sts.stats.stlo,eq_lat,eq_lon)[0]
for sts in stNs:
    sts.stats.distance = gps2dist_azimuth(sts.stats.stla,sts.stats.stlo,eq_lat,eq_lon)[0]
for sts in stZs:
    sts.stats.distance = gps2dist_azimuth(sts.stats.stla,sts.stats.stlo,eq_lat,eq_lon)[0]

# 绘制各方向的Record_Section图
stEs.plot(type='section', outfile = outputdir + '/EIE_Section.png')
stNs.plot(type='section', outfile = outputdir + '/EIN_Section.png')
stZs.plot(type='section', outfile = outputdir + '/EIZ_Section.png')

time2 = time.time()
print("各方向上的section图绘制完毕！花费时间为%d秒!"  % (time2 - time1))

# 对各方向的波形除以系数，转换为加速度文件
for ssts in stEs:
    ssts.data = ssts.data / ratio

for ssts in stNs:
    ssts.data = ssts.data / ratio

for ssts in stZs:
    ssts.data = ssts.data / ratio

# 把地震事件前10s的平均值作为矫正值
for sts in stEs:
    cor = np.mean(sts.data[0:1000])
    sts.data = sts.data - cor
for sts in stNs:
    cor = np.mean(sts.data[0:1000])
    sts.data = sts.data - cor
for sts in stZs:
    cor = np.mean(sts.data[0:1000])
    sts.data = sts.data - cor

## 计算理论到时，把理论到时之前10s的平均值作为矫正值
## iasp91模型 
#model = TauPyModel(model="iasp91")
#for sts in stEs:
#    arrivals = model.get_travel_times_geo(source_depth_in_km = eq_dep, source_latitude_in_deg = eq_lat, source_longitude_in_deg = eq_lon, receiver_latitude_in_deg = sts.stats.stla, receiver_longitude_in_deg = sts.stats.stlo)
#    arr = arrivals[0]
#    cor_st = int((arr.time + 50)*sts.stats.sampling_rate)
#    cor_en = int((arr.time + 60)*sts.stats.sampling_rate)
#    cor = np.mean(sts.data[cor_st:cor_en])
#    sts.data = sts.data - cor
#for sts in stNs:
#    arrivals = model.get_travel_times_geo(source_depth_in_km = eq_dep, source_latitude_in_deg = eq_lat, source_longitude_in_deg = eq_lon, receiver_latitude_in_deg = sts.stats.stla, receiver_longitude_in_deg = sts.stats.stlo)
#    arr = arrivals[0]
#    cor_st = int((arr.time + 50)*sts.stats.sampling_rate)
#    cor_en = int((arr.time + 60)*sts.stats.sampling_rate)
#    cor = np.mean(sts.data[cor_st:cor_en])
#    sts.data = sts.data - cor
#for sts in stZs:
#    arrivals = model.get_travel_times_geo(source_depth_in_km = eq_dep, source_latitude_in_deg = eq_lat, source_longitude_in_deg = eq_lon, receiver_latitude_in_deg = sts.stats.stla, receiver_longitude_in_deg = sts.stats.stlo)
#    arr = arrivals[0]
#    cor_st = int((arr.time + 50)*sts.stats.sampling_rate)
#    cor_en = int((arr.time + 60)*sts.stats.sampling_rate)
#    cor = np.mean(sts.data[cor_st:cor_en])
#    sts.data = sts.data - cor

time3 = time.time()
print("加速度波形校正完毕！花费时间为%d秒!"  % (time3 - time2))
    
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
    EW = sts_E[0].max()
    NS = sts_N[0].max()
    UD = sts_Z[0].max()
    thr = np.max(np.sqrt(np.square(sts_E[0].data) + np.square(sts_N[0].data) + np.square(sts_Z[0].data)))
    if thr == 0:
        i_pga = 0
    else:
        i_pga = 3.17*math.log10(thr) + 6.59
    pga_EW.append(EW)
    pga_NS.append(NS)
    pga_UD.append(UD)
    pga_thr.append(thr)
    I_pga.append(i_pga)

time4 = time.time()
print("PGA计算完毕！花费时间为%d秒!"  % (time4 - time3))
    
#以下开始进行pgv的计算
stEs_v = stEs.integrate()
stEs_v.filter('bandpass', freqmin=1.0, freqmax=10)
stNs_v = stNs.integrate()
stNs_v.filter('bandpass', freqmin=1.0, freqmax=10)
stZs_v = stZs.integrate()
stZs_v.filter('bandpass', freqmin=1.0, freqmax=10)
    
#按台站名分别计算各个方向的pgv和合成pgv
pgv_EW = []
pgv_NS = []
pgv_UD = []
pgv_thr = []
I_pgv = []
for sta in stas:
    sts_E = stEs_v.select(station = sta)
    sts_N = stNs_v.select(station = sta)
    sts_Z = stZs_v.select(station = sta)
    EW = sts_E[0].max()
    NS = sts_N[0].max()
    UD = sts_Z[0].max()
    thr = np.max(np.sqrt(np.square(sts_E[0].data) + np.square(sts_N[0].data) + np.square(sts_Z[0].data)))
    if thr == 0:
        i_pgv = 0
    else:
        i_pgv = 3.00*math.log10(thr) + 9.77
    pgv_EW.append(EW)
    pgv_NS.append(NS)
    pgv_UD.append(UD)
    pgv_thr.append(thr)
    I_pgv.append(i_pgv)
       
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

time5 = time.time()
print("PGV和仪器烈度计算完毕！花费时间为%d秒!"  % (time5 - time4))
        
#到此计算完毕，以下为结果输出 
#取出经纬度，震中距,台站名等信息
st_lat = []
st_lon = []
distance = []
st_loc = []
for sta in stas:
    st_lat.append(stEs.select(station = sta)[0].stats.stla)
    st_lon.append(stEs.select(station = sta)[0].stats.stlo)
    distance.append((stEs.select(station = sta)[0].stats.distance)/1000)   
    st_loc.append((stEs.select(station = sta)[0].stats.loc))
      
# 把结果写入到文件
r1 = np.array([stas, st_loc, st_lat, st_lon, distance, pga_EW, pga_NS, pga_UD, pga_thr, pgv_EW, pgv_NS, pgv_UD, pgv_thr, I])
r1 = np.ndarray.tolist(r1.T)
result = sorted(r1, key=lambda r1_tuple: float(r1_tuple[4]))
name = ['station', 'location', 'lat', 'lon', 'distance(km)', 'EW-PGA', 'NS-PGA', 'UD-PGA', 'PGA', 'EW-PGV', 'NS-PGV', 'UD-PGV', 'PGV', 'I']
intensity = pd.DataFrame(columns=name, data=result)
intensity.set_index(["station"], inplace=True)
intensity.to_csv(outputdir + '/py_intensity.csv')
# 另存为intensity.txt
with open(outputdir + '/py_intensity.csv','r', encoding='UTF-8') as r:
    lines=r.readlines()
with open(outputdir + '/intensity.txt','w',encoding="utf-8") as f_w:
    for line in lines:
        if "," in line:
            line = line.replace(",", " ")
        f_w.write(line)

##输出intensitydata和epicenterdata用于gmt作图
#with open(outputdir + '/py_intensity.csv','r') as r:
#    lines=r.readlines()
#with open('intensitydata','w') as w: 
#    w.writelines(lines[1:])
#    
#with open('intensitydata','r') as r:
#    lines=r.readlines()
#with open('intensitydata','w',encoding="utf-8") as f_w:
#    for line in lines:
#        if "," in line:
#            line = line.replace(",", " ")
#        f_w.write(line)
#        
#with open('epicenterdata','w') as w:
#    w.write(str(eq_year) + '/'+  str(eq_month) + '/'+  str(eq_date) + ' ' + 'M' + str(eq_mag) +' ' + str(eq_lat) + ' ' + str(eq_lon) + ' ' + str(eq_dep))
#
##调用gmt画图
#os.system('./pga.gmt intensitydata epicenterdata')
#os.system('./pgv.gmt intensitydata epicenterdata')
#os.system('./intensity.gmt intensitydata epicenterdata')
#os.system('ps2raster -A -P pga.ps')
#os.system('ps2raster -A -P pgv.ps')
#os.system('ps2raster -A -P intensity.ps')
#shutil.move('pga.jpg',outputdir + '/pga.jpg')
#shutil.move('pgv.jpg',outputdir + '/pgv.jpg')
#shutil.move('intensity.jpg',outputdir + '/intensity.jpg')
#shutil.move('pga.ps',outputdir + '/pga.ps')
#shutil.move('pgv.ps',outputdir + '/pgv.ps')
#shutil.move('intensity.ps',outputdir + '/intensity.ps')
#os.remove('intensitydata')
#os.remove('epicenterdata')

# 获取结束时间
time_over = time.time()
print("计算完毕！总共花费时间为%d秒!"  % (time_over - time_begin))