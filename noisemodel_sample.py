# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:46:03 2020

@author: Kaiyi Jiang
"""
import pandas as pd
import numpy as np


def binary_search(arr, low, high, x):
    max_value = arr.max()
    min_value = arr.min()
    if x >= max_value:
        return len(arr) - 2
    elif x < min_value:
        return 0
    else:
        mid = round((high + low) / 2)

        if x >= arr[mid] and x < arr[mid + 1]:
            return mid

        elif x < arr[mid]:
            return binary_search(arr, low, mid - 1, x)
        elif x >= arr[mid + 1]:
            return binary_search(arr, mid + 1, high, x)


def noiselevel(stain, GFP, xbins, ybins):
    # calculate probablity of real given observed
    # Data stored with stain being observed, GFP being real
    combine = pd.concat([stain, GFP], axis=1)
    combine = np.log10(combine)
    x = combine.iloc[:, 0]
    y = combine.iloc[:, 1]
    xmin = x.min()
    ymin = y.min()
    xmax = x.max()
    ymax = y.max()
    xbin = np.zeros(xbins + 1)
    ybin = np.zeros(ybins + 1)
    xbin[0] = xmin
    xbin[xbins] = xmax
    ybin[0] = ymin
    ybin[ybins] = ymax

    for index in range(xbins - 1):
        diff = (xmax - xmin) / (xbins)
        add = xbin[index] + diff
        xbin[index + 1] = add

    for index in range(ybins - 1):
        diff = (ymax - ymin) / (ybins)
        add = ybin[index] + diff
        ybin[index + 1] = add

    histogram = np.zeros((ybins, xbins))

    for index in range(xbins):
        if index < xbins - 1:
            mask = combine.iloc[:, 0] >= xbin[index]
            mask2 = combine.iloc[:, 0] < xbin[index + 1]
        else:
            mask = combine.iloc[:, 0] >= xbin[index]
            mask2 = combine.iloc[:, 0] <= xbin[index + 1]

        gfp = combine[mask & mask2]
        gfp = gfp.iloc[:, 1]
        count, bin_edge = np.histogram(gfp, bins=ybin)
        summation = sum(count)
        if summation == 0:
            histogram[:, index] = histogram[:, index - 1]
        else:
            count = count / summation
            histogram[:, index] = count
    xbin = 10 ** (xbin)
    ybin = 10 ** (ybin)
    list = [histogram, xbins, xbin, ybins, ybin]
    return list


def findreal(list, kinase):
    histogram = list[0]
    xbins = list[1]
    xbin = list[2]
    ybins = list[3]
    ybin = list[4]
    low = 0
    high = xbins - 1

    k_index = binary_search(xbin, low, high, kinase)

    pdf = histogram[:, k_index]

    mean_y = []
    for index in range(ybins):
        mean = (ybin[index] + ybin[index + 1]) / 2
        mean_y.append(mean)

    return pdf, mean_y

def sample(noise,item):
    pdf,mean_y=findreal(noise,item)
    summation=sum(pdf)
    if summation==0:
        draw=0
    else:
        draw=np.random.choice(mean_y,100,True,pdf)
    return draw

def singlecell(list_k, list_s, item):
    kinase=item[0]
    substrate=item[1]
    PE=item[2]
    factor=item[3:len(item)]


    k = sample(list_k,kinase)
    s = sample(list_s,substrate)


    combine=np.vstack([k,s])
    combine=combine.T

    merge=np.tile(factor,(len(combine),1))
    combine=np.hstack([combine,merge])
    return combine


# matrix_noise = pd.read_csv('FP-tag data/PEnoise.csv')
#
# matrix_noise = matrix_noise[matrix_noise['hospho : PE - Area'] > 0]
# matrix_noise = matrix_noise[matrix_noise['mCherry - Area'] > 0]
# pe = matrix_noise['hospho : PE - Area']
# mcherry = matrix_noise['mCherry - Area']
# xbins=100
# ybins=100
# list_pe=noiselevel(pe, mcherry, xbins, ybins)
#
# matrix_noise = pd.read_csv('FP-tag data/Flagnoise.csv')
#
# matrix_noise = matrix_noise[matrix_noise['Flag : Alexa Fluor 680 - Area'] > 0]
# matrix_noise = matrix_noise[matrix_noise['GFP - Area'] > 0]
# flag = matrix_noise['Flag : Alexa Fluor 680 - Area']
# gfp = matrix_noise['GFP - Area']
# list_k=noiselevel(flag, gfp, xbins, ybins)
#
# matrix_noise = pd.read_csv('FP-tag data/Mycnoise.csv')
#
# matrix_noise = matrix_noise[matrix_noise['Myc : Alexa Fluor 750 - Area'] > 0]
# matrix_noise = matrix_noise[matrix_noise['GFP - Area'] > 0]
# myc = matrix_noise['Myc : Alexa Fluor 750 - Area']
# gfp = matrix_noise['GFP - Area']
# list_s=noiselevel(myc, gfp, xbins, ybins)
#
#
#
# data = pd.read_csv('Globalfit_data/170_127.csv')
#
# data = data[data['Flag : Alexa Fluor 680 - Area'] > 0]
# data = data[data['Myc : Alexa Fluor 750 - Area'] > 0]
# data = data[data['hospho : PE - Area'] > 0]
# data1 = data
# kinase = data['Flag : Alexa Fluor 680 - Area']
# substrate = data['Myc : Alexa Fluor 750 - Area']
# PE = data['hospho : PE - Area']
# item=np.arange(4)
# a=singlecell(list_k, list_s, 100,1,0.2, item)





