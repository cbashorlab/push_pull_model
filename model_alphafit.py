# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:46:03 2020

@author: Kaiyi Jiang
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import concurrent.futures 
import time

tic=time.perf_counter()
matrix_noise=pd.read_csv('Staining_noise.csv')
data=pd.read_csv('Kinase_Sub_P.csv')

data=data[data['Flag : Alexa Fluor 680 - Area']>0]
data=data[data['Myc : Alexa Fluor 750 - Area']>0]
data=data[data['hospho : PE - Area']>0]
data1=data
data=data.head(10)
kinase=data['Flag : Alexa Fluor 680 - Area']
substrate=data['Myc : Alexa Fluor 750 - Area']
PE=data['hospho : PE - Area']


matrix_noise=matrix_noise[matrix_noise['CD8 : APC - Area']>0]
matrix_noise=matrix_noise[matrix_noise['GFP - Area']>0]
myc=matrix_noise['CD8 : APC - Area']
GFP=matrix_noise['GFP - Area']

def binary_search(arr,low,high, x):
    max_value = arr.max()
    min_value=arr.min()
    if x>=max_value:
        return len(arr)-2
    elif x<min_value:
        return 0
    else:
        mid = round((high + low) / 2)

        if x>=arr[mid] and x<arr[mid+1]:
            return mid

        elif x<arr[mid]:
            return binary_search(arr, low, mid - 1, x)
        elif x>=arr[mid+1]:
            return binary_search(arr, mid + 1, high, x)

def noiselevel(stain,GFP,xbins,ybins):
    #calculate probablity of real given observed
    #Data stored with stain being observed, GFP being real
    combine=pd.concat([stain,GFP],axis=1)
    combine=np.log10(combine)
    x=combine.iloc[:,0]
    y=combine.iloc[:,1]
    xmin=x.min()
    ymin=y.min()
    xmax=x.max()
    ymax=y.max()
    xbin=np.zeros(xbins+1)
    ybin=np.zeros(ybins+1)
    xbin[0]=xmin
    xbin[xbins]=xmax
    ybin[0]=ymin
    ybin[ybins]=ymax
    
    for index in range(xbins-1):
        diff=(xmax-xmin)/(xbins)
        add=xbin[index]+diff
        xbin[index+1]=add
    
    for index in range(ybins-1):
        diff=(ymax-ymin)/(ybins)
        add=ybin[index]+diff
        ybin[index+1]=add
    
    histogram=np.zeros((ybins,xbins))
    
    for index in range(xbins):
        if index < xbins-1:
            mask=combine.iloc[:,0]>=xbin[index]
            mask2=combine.iloc[:,0]<xbin[index+1]
        else:
            mask=combine.iloc[:,0]>=xbin[index]
            mask2=combine.iloc[:,0]<=xbin[index+1]

        gfp=combine[mask & mask2]
        gfp=gfp.iloc[:,1]
        count,bin_edge=np.histogram(gfp,bins=ybin)
        summation=sum(count)
        count=count/summation
        histogram[:,index]=count
    xbin=10**(xbin)
    ybin=10**(ybin)
    list=[histogram,xbins,xbin,ybins,ybin]
    return list

noise_phosphotase=noiselevel(myc,GFP,100,100)


def findreal(list,kinase):
    histogram=list[0]
    xbins=list[1]
    xbin=list[2]
    ybins=list[3]
    ybin=list[4]
    low=0
    high=xbins-1

    k_index=binary_search(xbin,low,high,kinase)

    pdf=histogram[:,k_index]
    
    mean_y=[]
    for index in range(ybins):
        mean=(ybin[index]+ybin[index+1])/2
        mean_y.append(mean)
    
    return pdf,mean_y


    
    
def singlecell(list_k,list_s,list_pe,kinase,substrate,PE):
    k_pdf,mean_k=findreal(list_k,kinase)
    s_pdf,mean_s=findreal(list_s,substrate)
    p_pdf,mean_p=findreal(list_pe,PE)


    prob=np.stack(np.meshgrid(k_pdf,s_pdf,p_pdf),-1).reshape(-1,3)
    alpha_store=np.stack(np.meshgrid(mean_k,mean_s,mean_p),-1).reshape(-1,3)
    k=alpha_store[:,0]
    s=alpha_store[:,1]
    p=alpha_store[:,2]
    probablity=prob[:,0]*prob[:,1]*prob[:,2]
    alpha_calc=(p-k)*(p-s)/p

    combine=np.vstack([alpha_calc,probablity])
    return combine


def multicell(list_k,list_s,list_pe,kinase,substrate,PE):

    combine=pd.concat([kinase,substrate,PE],axis=1)
    size=len(combine)
    total = np.zeros((2, list_k[3] ** 3 * size))
    for index in range(size):
        a=singlecell(list_k,list_s,list_pe,combine.iloc[index,0],combine.iloc[index,1],combine.iloc[index,2])
        start=list_k[3]**3*index
        end=list_k[3]**3*(index+1)
        total[:,start:end]=a

    # def helper(kinase,substrate,PE):
    #     return singlecell(histogram_k,histogram_s,histogram_p,xbins,xbin,ybins,ybin,kinase,substrate,PE)
    #
    # total=np.array(([0],[0]))
    # if __name__ == '__main__':
    #     with concurrent.futures.ThreadPoolExecutor() as executors:
    #       for result in executors.map(helper,kinase,substrate,PE):
    #           total=np.hstack([total,result])
    # total=np.delete(total,0,axis=1)

    return total

def alphacalc(alpha_bins,list_k,list_s,list_pe,kinase,substrate,PE):
              
    total=multicell(list_k,list_s,list_pe,kinase,substrate,PE)
    #how to process alpha 
    #for now alpha must be positive
    # log probablity is used
    total=total.T
    total=total[total[:,1]>0]
    total=total[total[:,0]>0]
    total=total.T
    total[1,:]=np.log10(total[1,:])
    
    min_1=total[0,:].min()
    max_1=total[0,:].max()
    min_1=np.log10(min_1)
    max_1=np.log10(max_1)
    bin_alpha=np.zeros(alpha_bins+1)
    bin_alpha[0]=min_1
    bin_alpha[alpha_bins]=max_1
    for index in range(alpha_bins-1):
        diff=(max_1-min_1)/(alpha_bins)
        add=bin_alpha[index]+diff
        bin_alpha[index+1]=add
    bin_alpha=10**(bin_alpha)
    
    mean_alpha=[]
    for index in range(alpha_bins):
        mean=(bin_alpha[index]+bin_alpha[index+1])/2
        mean_alpha.append(mean)
    prob_alpha=np.zeros(len(mean_alpha))


    for n in range(len(total.T)):
        index=binary_search(bin_alpha,0,alpha_bins-1,total[0,n])
        prob_alpha[index]=prob_alpha[index]+total[1,n]
    
    prob_alpha=prob_alpha/sum(prob_alpha)
    
    return mean_alpha, prob_alpha

mean_alpha, prob_alpha=alphacalc(20,noise_phosphotase,noise_phosphotase,noise_phosphotase,kinase,substrate,PE)


plt.plot(mean_alpha,prob_alpha)
plt.xlabel('Alpha')
plt.ylabel('Normalized Probablity')
plt.title('Alpha Probablity')
plt.xscale('log')
plt.show()
max_index=np.argmax(prob_alpha)
alpha_max=mean_alpha[max_index]
print(alpha_max)
toc=time.perf_counter()
print(toc-tic)


def simulation(mean_alpha,prob_alpha,data):
    kinase=data['Flag : Alexa Fluor 680 - Area']
    #kinase=np.log10(kinase)
    substrate=data['Myc : Alexa Fluor 750 - Area']
    #substrate=np.log10(substrate)
    PE=data['hospho : PE - Area']
    #PE=np.log10(PE)
    max_index=np.argmax(prob_alpha)
    alpha_max=mean_alpha[max_index]
    print(alpha_max)
    pe_predicted=[]
    for index,value in PE.items():
       sum = kinase[index]+ substrate[index] + alpha_max
       if ((sum**2)-4*substrate[index]*kinase[index]) <=0:
           a=0.5*sum
       else:
           a=0.5*(sum+((sum**2)-4*substrate[index]*kinase[index])**0.5)
           if a>substrate[index]:
               a=0.5*(sum-((sum**2)-4*substrate[index]*kinase[index])**0.5)
       pe_predicted.append(a)
    return pe_predicted

# pe_predicted=simulation(mean_alpha,prob_alpha,data1)
# fig,(ax1,ax2)=plt.subplots(1,2)
# ax1.scatter(np.log10(data1['Flag : Alexa Fluor 680 - Area']),np.log10(data1['hospho : PE - Area']))
# ax2.scatter(np.log10(data1['Flag : Alexa Fluor 680 - Area']),np.log10(pe_predicted))
# ax1.set_title('Original')
# ax2.set_title('Predicted')
# ax1.set_xlabel('Kinase')
# ax1.set_ylabel('Phosphorylation')
# ax2.set_xlabel('Kinase')

# plt.show()
