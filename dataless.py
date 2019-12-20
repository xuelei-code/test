# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 15:48:47 2018

@author: lei
"""
import obspy
from obspy import read, read_inventory
from os import listdir
from obspy.core import UTCDateTime


# 读入时间
t0 = UTCDateTime("2018-04-22T16:30:00.000")
# 读入仪器响应文件
inv = read_inventory('F:/work/mycode/test/FJ.dataless.xml')
pre_filt = [0.001, 0.1, 10, 20]
# 读入波形数据
filepath='F:/work/mycode/test/mseed'
filename_list=listdir(filepath)

for filename in filename_list:
    st = read(filepath + '/' + filename,starttime = t0 - 60, endtime = t0 + 300)

    # the corresponding response is included in ObsPy as a StationXML file

    # the routine automatically picks the correct response for each trace
    # define a filter band to prevent amplifying noise during the deconvolution
    re_filt = (0.005, 0.1, 10, 11)
    Vst = st.remove_response(inventory=inv, output='VEL', pre_filt=pre_filt)
    Ast = st.remove_response(inventory=inv, output='ACC', pre_filt=pre_filt)
    # Vst.plot()
    # Ast.plot()
    tr = Vst[0]
    print(tr.stats)