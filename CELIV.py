# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 10:41:17 2021
Pick the peak od CELIV from signal data extracted from background
Program work best when the square impulse in dark is at the beginning of the half of the signal
program will get error when the square impulse in dark is far end of the timescale
squarepulse must be at least 50 data points
Tested in Python3.9
@author: Kai Lin WOON
"""
import csv
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import os

def reading_file(path):    
    file=open(path)
    csvreader = csv.reader(file)
    x,y = [],[]
    for row in csvreader:
        x0=float(row[0])
        y0=float(row[1])
        x.append(x0)
        y.append(y0)
    file.close()
    return x,y

def finding_onset_fall(x,y):
    def take_care_of_shift_peak(x,y,n):
        n=int(len(y)/2)+n  # look at the beginnging half
        myindex=[]
        findingonset=[]
        for i in range(n):
            findingonset.append(abs(y[i+1]-y[i]))
        cutoff=max(findingonset)
        for i in range(n):
            if abs(y[i+2]-y[i])>0.5*cutoff:
                myindex.append(i)
        t_index0=myindex[0]+1
        myindex=[]
        for i in range(t_index0+15,len(y)-2):
            if abs(y[i+2]-y[i])>0.5*cutoff:
                myindex.append(i)
        t_index1=myindex[0]-1
        return t_index0, t_index1
    n=0
    t_index0,t_index1=take_care_of_shift_peak(x,y,n)
    while t_index1-t_index0<50:  #if the interval is too short, look at later laft of the signal
        n=n+50
        t_index0, t_index1=take_care_of_shift_peak(x,y,n)
    print( t_index0, t_index1)
    return t_index0,t_index1

    
def format_data(x,y,t_index0,t_index1):    
    newx,newy=[],[]
    for i in range(len(x)):
      if t_index1>i>t_index0:
          newx.append(x[i])
          newy.append(y[i])
    return newx,newy
          
def substract_dark(x0,y0,y1,filter_n):
    x0=np.array(x0)
    y0=np.array(y0)
    y1=np.array(y1)
    y_unfiltered=y1-y0
    y_filtered = savgol_filter(y_unfiltered, filter_n+21, 3) #smoothing parameters
    return x0,y_filtered,y_unfiltered



def finding_transit_time(path1,path0):
    x1,y1=reading_file(path1)
    x0,y0=reading_file(path0)
    t_index0,t_index1=finding_onset_fall(x0,y0)  # use dark
    newx1,newy1=format_data(x1,y1,t_index0,t_index1)  #cut necessary data
    newx0,newy0=format_data(x0,y0,t_index0,t_index1)  #cut necessary data
    
    mydata=[]
    for i in range(14):
        newx0,y_filtered,y_unfiltered=substract_dark(newx0,newy0,newy1,i*4)
        location=np.argmax(y_filtered)
        transit_time=newx0[location]-newx0[0]
        mydata.append(transit_time)
    return np.median(mydata), t_index0,t_index1

directory='C:/Users/user/Downloads/CELIV_new_10.09.2021/CELIV_new_10.09.2021' #where is CELIV file only CSV file
arr = os.listdir(directory) 
f=open('C:/Users/user/Downloads/test2.txt','w')  #write file

mydata=[]
for i in range(int(len(arr)/2)):
    file_even=arr[i*2]
    file_odd=arr[i*2+1]
    path0=directory+'/'+file_even
    path1=directory+'/'+file_odd
    print(file_even)
    try:
        transit,t_index0,t_index1=finding_transit_time(path1,path0)
    except:
        transit='error'
    f.write(str(file_even)+' , '+str(transit)+' , '+str(t_index0)+' , '+str(t_index1))
    f.write('\n')
f.close()

# fig, ax = plt.subplots()
# ax.scatter(newx0,y_unfiltered)
# ax.plot(newx0,y_filtered, 'r')
# ax.set(xlim=(0,np.max(newx1)),ylim=(0,np.max(y_filtered)))
# plt.show()     
