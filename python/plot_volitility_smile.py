# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:30:07 2016

@author: mahone
"""
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt
import numpy as np
import xlrd

workbook = xlrd.open_workbook('spxopt.xlsx')
first_sheet = workbook.sheet_by_index(0)
m =[float( first_sheet.cell_value(i, 2)) for i in range(1,first_sheet.nrows)]
t =[float( first_sheet.cell_value(i, 1))for i in range(1,first_sheet.nrows)]
k =[float( first_sheet.cell_value(i, 3))for i in range(1,first_sheet.nrows)]
v =[first_sheet.cell_value(i, 5)for i in range(1,first_sheet.nrows)]
v = np.asarray(v)
m = np.asarray(m)
k = np.asarray(k)
v[v=='']='0'
v = v.astype(np.float)
for i in range(586):
            m[i] = m[i] - t[i]

x, y = np.meshgrid(k, m)
fig = plt.figure(figsize=(20, 10))
plot = p3.Axes3D(fig)
plot.plot_wireframe(x, y, v)
plot.set_xlabel('Strike Price')
plot.set_ylabel('TTM')
plot.set_zlabel('Impled volatility')