import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import shgo
import matplotlib.pyplot as plt

def onehotencoding(list,list2,number):
    vec = pd.get_dummies(list)
    vec = vec.values
    list2.append('A')
    vec2 = pd.get_dummies(list2)
    vec2 = vec2.values
    del list2[len(list2) - 1]

    default = list[0]
    data = pd.read_csv(default)
    data = data.sample(number)
    data = data.reset_index(drop=True)
    factor_zipper = pd.Series([default] * len(data), name='vk')
    factor_velocity = pd.Series(['A'] * len(data), name='vp')
    data = pd.concat([data, factor_zipper, factor_velocity], axis=1)

    del list[0]
    for name in list:
        data1 = pd.read_csv(name)
        data1 = data1.sample(number)
        data1 = data1.reset_index(drop=True)
        factor_zipper = pd.Series([name] * len(data1), name='vk')
        factor_velocity = pd.Series(['A'] * len(data1), name='vp')
        data1 = pd.concat([data1, factor_zipper,factor_velocity], axis=1)
        data = pd.concat([data, data1])

    name = list2[0]
    data1 = pd.read_csv(name)
    data1 = data1.sample(number)
    data1 = data1.reset_index(drop=True)
    factor_zipper = pd.Series([default] * len(data1), name='vk')
    factor_velocity = pd.Series([name] * len(data1), name='vp')
    data1 = pd.concat([data1, factor_zipper, factor_velocity], axis=1)
    data = pd.concat([data, data1])

    del list2[0]
    for name in list2:
        data1 = pd.read_csv(name)
        data1 = data1.sample(number)
        data1 = data1.reset_index(drop=True)
        factor_zipper = pd.Series([default] * len(data1), name='vk')
        factor_velocity = pd.Series([name] * len(data1), name='vp')
        data1 = pd.concat([data1, factor_zipper,factor_velocity], axis=1)
        data = pd.concat([data, data1])


    factor = pd.get_dummies(data.vk)
    factor1 = pd.get_dummies(data.vp)
    data = data.drop(columns=['vk'])
    data = data.drop(columns=['vp'])
    data = pd.concat([data, factor, factor1], axis=1)

    return data,vec,vec2

list=['k10p10.csv','k5p10.csv','k15p10.csv']
list2=['k10p5.csv','k10p15.csv']
data,vec,vec_p=onehotencoding(list,list2,6000)
list=['k10p10.csv','k5p10.csv','k15p10.csv']
list2=['k10p5.csv','k10p15.csv']


kinase=data['Kinase']
substrate=data['Substrate']
pe=data['PE']
phosphatase=data['Phosphatase']
factor=data.iloc[:,5:data.shape[1]]

kinase=kinase.values
substrate=substrate.values
pe=pe.values
phosphatase=phosphatase.values
factor=factor.values

vec1=vec[0,:]
vec2=vec[1,:]
vec3=vec[2,:]

vec4=vec_p[2,:]
vec5=vec_p[0,:]
vec6=vec_p[1,:]

data=np.vstack([kinase,substrate,phosphatase,pe,factor.T])



def costfunc(input):
    alpha,alpha2,v_k1,v_k2,v_k3,v_p1,v_p2,v_p3=input
    vbg=0
    k=data[0,:]
    s=data[1,:]
    p=data[2,:]
    pe=data[3,:]
    i = data[4:(len(list) + 4), :]
    j = data[(len(list) + 4):data.shape[0], :]

    kinase=np.ones(k.shape)
    phosphatase=np.ones(k.shape)
    old_substrate=np.ones(k.shape)
    substrate=s/(kinase/alpha+phosphatase/alpha2+1)
    loop=0
    while np.sum(abs((substrate-old_substrate)/old_substrate))>0.01:
        loop=loop+1
        old_substrate=substrate
        kinase=k/(1+(substrate/alpha))
        phosphatase=p/(1+(substrate/alpha2))
        substrate = s / (kinase / alpha + phosphatase / alpha2 + 1)
    predict=s*(((np.dot(vec1,i)*v_k1)+(np.dot(vec2,i)*v_k2)+(np.dot(vec3,i)*v_k3))*kinase/alpha + vbg)/(((np.dot(vec1,i)*v_k1)+(np.dot(vec2,i)*v_k2)+(np.dot(vec3,i)*v_k3))*kinase/alpha + vbg + 1 + ((np.dot(vec4,j)*v_p1)+(np.dot(vec5,j)*v_p2)+(np.dot(vec6,j)*v_p3))*phosphatase/alpha2)
    cost=predict-pe
    cost=np.power(cost,2)
    cost=np.sum(cost)
    print(loop)
    #print(np.sum(np.power(cost,2)))
    return cost
# bound1=((0,None),(0,None),(0,None),(0,None))
res=minimize(costfunc,[1000,1000,10,10,10,10,10,10],method='BFGS')
print(res)
# bound1=([0,0,0,0],[10**15,10**15,10**15,10**15])
# bound2=((0.1,None),(0.1,None),(0.1,None),(0.1,None))
# res=minimize(costfunc,[1000,1000,100,10])
# print(res)