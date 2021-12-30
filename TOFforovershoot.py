"""
Created on 26/12/2021 2021
@author: Dr KL Woon
Time of Flight
Finding transit time by finding the maximum angle
between two lines (here negative angle)
the lines are best-fitted linear line in log-log curve
Only R2 above certain values are used. Interval of data to find
the best fitted lines mimics the normal way doing
Note very long oscillitory overshots might give incorrect result
The onset of photocurrrent must be at zero time
Tested in Python3.9
"""
import csv
import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import os     

def fitfunction(func,x,y):
    popt, pcov = curve_fit(func, x, y)                                                             
    predicted=func(x, *popt) 
    r=np.corrcoef(y, predicted)
    r2=r[0][1]**2
    return popt, r2

def func(m,x,c):
    return m*x+c

def reading_file(path):    
    file=open(path)
    csvreader = csv.reader(file)
    x,y = [],[]
    for row in csvreader:
        x0=float(row[0])
        y0=float(row[1])
        if x0>=0:    # data saved when x is greater than zero
            x.append(x0)
            y.append(y0)
    file.close()
    x=np.array(x)
    y=np.array(y)
    y=(y-min(y))*max(x)/max(y)  #make x and y are equal propotion
    y_filtered = savgol_filter(y, 17, 3) #smoothing parameters
    y_filtered[y_filtered<0] = 1.0e-09
    return x,y_filtered

def createinterval(totaldatapoint,interval):
    x0=totaldatapoint
    x=totaldatapoint
    n,i=1,0
    mylist=[]
    while x>0:
        x=x-interval*n
        mylist.append(interval*n)
        i=i+1
        if i%2==0 and i!=0:
            n=n+1
    mylist=mylist[:-1]
    mylist.sort(reverse=True)
    mylist[0]=mylist[0]+(x0-sum(mylist))
    mylist2=[mylist[0]]
    for i in range(len(mylist)-1):
        z=sum(mylist[0:i+2])
        mylist2.append(z)
    return mylist2

def search_transit(mylist2,y_filtered,x):
    gradient,c,ang,rr2=[],[],[],[]
    ii=0
    for i in range(len(mylist2)-1):  #start from last points and progressing to beginning
        if ii==0:
            newy=y_filtered[-mylist2[i]:]
            newx=x[-1*mylist2[i]:]
        if ii!=0:
            newy=y_filtered[-mylist2[i+1]:-mylist2[i]]
            newx=x[-mylist2[i+1]:-mylist2[i]]
        ii=ii+1
        newy[newy==0]=0.0000001
        newy=np.log10(newy)
        newx[newx==0]=0.0000001
        newx=np.log10(newx)
        popt, r2=fitfunction(func,newx,newy)
        rr2.append(r2)
    ii=0
    for i in range(len(mylist2)-1):  #start from last points and progressing to beginning
        if ii==0:
            newy=y_filtered[-mylist2[i]:]
            newx=x[-1*mylist2[i]:]
        if ii!=0:
            newy=y_filtered[-mylist2[i+1]:-mylist2[i]]
            newx=x[-mylist2[i+1]:-mylist2[i]]
        ii=ii+1
        newy=np.log10(newy)
        newx[newx==0]=0.0000001
        newx=np.log10(newx)
        popt, r2=fitfunction(func,newx,newy)
        cut_off=np.median(rr2)*0.9
        if r2>cut_off:
            gradient.append(popt[0])
            c.append(popt[1])
    #scale not the same
    for i in range(len(c)-1):
        angle=np.arctan((gradient[i]-gradient[i+1])/(1+gradient[i+1]*gradient[i])) #calculate angle between 2 lines
        ang.append(angle*180/np.pi)
        # print(angle*180/np.pi)
    
    min_angle=min(ang) # assumption when angle is positive
    location=ang.index(min(ang))     
    m1=gradient[location]
    m2=gradient[location+1]
    c1=c[location]
    c2=c[location+1]
    intersect_x=(c2-c1)/(m1-m2)
    return np.power(10,intersect_x)

directory='C:/Users/user/Downloads/HD73_02.09.2021/HD73_02.09.2021' #where is TOF file only CSV file
arr = os.listdir(directory) 
f=open('C:/Users/user/Downloads/test.txt','w')  #write file
for file in arr:
    path=directory+'/'+file
    x,y_filtered=reading_file(path)
    mydata=[]
    for i in range(10):
        mylist2=createinterval(len(x),10+i) 
        intersect_x=search_transit(mylist2,y_filtered,x)
        mydata.append(intersect_x)
    to_write=np.median(mydata)
    print(to_write)
    f.write(str(file)+' , '+str(to_write))
    f.write('\n')
f.close()
