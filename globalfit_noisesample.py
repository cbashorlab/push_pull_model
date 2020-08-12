import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import noisemodel_sample
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import basinhopping
import time


#
# matrix_noise = pd.read_csv('FP-tag data/PEnoise.csv')
#
# matrix_noise = matrix_noise[matrix_noise['hospho : PE - Area'] > 0]
# matrix_noise = matrix_noise[matrix_noise['mCherry - Area'] > 0]
# pe = matrix_noise['hospho : PE - Area']
# mcherry = matrix_noise['mCherry - Area']

# list_pe=noisemodel_sample.noiselevel(pe, mcherry, xbins, ybins)
xbins=500
ybins=500
matrix_noise = pd.read_csv('FP-tag data/Flagnoise.csv')

matrix_noise = matrix_noise[matrix_noise['Flag : Alexa Fluor 680 - Area'] > 0]
matrix_noise = matrix_noise[matrix_noise['GFP - Area'] > 0]
flag = matrix_noise['Flag : Alexa Fluor 680 - Area']
gfp = matrix_noise['GFP - Area']
list_k=noisemodel_sample.noiselevel(flag, gfp, xbins, ybins)

matrix_noise = pd.read_csv('FP-tag data/Mycnoise.csv')

matrix_noise = matrix_noise[matrix_noise['Myc : Alexa Fluor 750 - Area'] > 0]
matrix_noise = matrix_noise[matrix_noise['GFP - Area'] > 0]
myc = matrix_noise['Myc : Alexa Fluor 750 - Area']
gfp = matrix_noise['GFP - Area']
list_s=noisemodel_sample.noiselevel(myc, gfp, xbins, ybins)

def onehotencoding(list,list2,number):
    vec = pd.get_dummies(list)
    vec = vec.values
    list2.append('A')
    vec2 = pd.get_dummies(list2)
    vec2 = vec2.values
    del list2[len(list2) - 1]

    default = list[0]
    data = pd.read_csv(default)
    data = data[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]
    data = data.sample(number)
    data = data.reset_index(drop=True)
    factor_zipper = pd.Series([default] * len(data), name='Zipper')
    factor_velocity = pd.Series(['A'] * len(data), name='Velocity')
    data = pd.concat([data, factor_zipper, factor_velocity], axis=1)

    del list[0]
    for name in list:
        data1 = pd.read_csv(name)
        data1 = data1[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]
        data1 = data1.sample(number)
        data1 = data1.reset_index(drop=True)
        factor_zipper = pd.Series([name] * len(data1), name='Zipper')
        factor_velocity = pd.Series(['A'] * len(data1), name='Velocity')
        data1 = pd.concat([data1, factor_zipper,factor_velocity], axis=1)
        data = pd.concat([data, data1])

    name = list2[0]
    data1 = pd.read_csv(name)
    data1 = data1[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]
    data1 = data1.sample(number)
    data1 = data1.reset_index(drop=True)
    factor_zipper = pd.Series([default] * len(data1), name='Zipper')
    factor_velocity = pd.Series([name] * len(data1), name='Velocity')
    data1 = pd.concat([data1, factor_zipper, factor_velocity], axis=1)
    data = pd.concat([data, data1])

    del list2[0]
    for name in list2:
        data1 = pd.read_csv(name)
        data1 = data1[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]
        data1 = data1.sample(number)
        data1 = data1.reset_index(drop=True)
        factor_zipper = pd.Series([default] * len(data1), name='Zipper')
        factor_velocity = pd.Series([name] * len(data1), name='Velocity')
        data1 = pd.concat([data1, factor_zipper,factor_velocity], axis=1)
        data = pd.concat([data, data1])

    data = data[data['Flag : Alexa Fluor 680 - Area'] > 0]
    data = data[data['Myc : Alexa Fluor 750 - Area'] > 0]
    data = data[data['hospho : PE - Area'] > 0]



    factor = pd.get_dummies(data.Zipper)
    factor1 = pd.get_dummies(data.Velocity)
    data = data.drop(columns=['Zipper'])
    data = data.drop(columns=['Velocity'])
    data = pd.concat([data, factor, factor1], axis=1)

    return data,vec,vec2

# List1 should include zipper files (the first file is default configuration)
# List 2 include kinase variants files
list1=['Globalfit_data/170_127.csv','Globalfit_data/170EEV.csv','Globalfit_data/170EES.csv','Globalfit_data/170EEK.csv','Globalfit_data/170RR.csv']
list2=['Globalfit_data/170L1.csv','Globalfit_data/170L2.csv']
data,vec,vec_va=onehotencoding(list1,list2,1000)
list1=['Globalfit_data/170_127.csv','Globalfit_data/170EEV.csv','Globalfit_data/170EES.csv','Globalfit_data/170EEK.csv','Globalfit_data/170RR.csv']
list2=['Globalfit_data/170L1.csv','Globalfit_data/170L2.csv']


kinase=data['Flag : Alexa Fluor 680 - Area']
substrate=data['Myc : Alexa Fluor 750 - Area']
PE=data['hospho : PE - Area']
factor=data.iloc[:,3:data.shape[1]]
kinase=kinase.values
substrate=substrate.values
PE=PE.values
factor=factor.values

#one hot for zipper
vec1=vec[0,:]
vec2=vec[1,:]
vec3=vec[2,:]
vec4=vec[3,:]
vec5=vec[4,:]
#one hot for va
vec6=vec_va[2,:]
vec7=vec_va[0,:]
vec8=vec_va[1,:]

tic=time.perf_counter()
combinedata=np.vstack([kinase,substrate,PE,factor.T])
noisedata=noisemodel_sample.singlecell(list_k, list_s, combinedata[:,0])
for index in range(combinedata.shape[1]-1):
    matrix=noisemodel_sample.singlecell(list_k, list_s, combinedata[:,index+1])
    noisedata=np.vstack([noisedata,matrix])
#Format: Kinase, Substrate, one hot factor
data=noisedata.T
ydata=np.zeros(len(data.T))
toc=time.perf_counter()
print(toc-tic)



def substratebound(data,alpha1,alpha2,alpha3,alpha4,alpha5,v1,v2,v3):
    vbg=0
    k=data[0,:]
    s=data[1,:]
    p=data[2,:]
    i=data[3:(len(list1)+3),:]
    j=data[(len(list1)+3):data.shape[0],:]


    sum = k + s + (np.dot(vec1,i)*alpha1) + (np.dot(vec2,i)*alpha2) + (np.dot(vec3,i)*alpha3) + (np.dot(vec4,i)*alpha4) + (np.dot(vec5,i)*alpha5)
    bound=0.5 * (sum - ((sum ** 2) - 4 * s * k) ** 0.5)
    #print(bound.shape)
    unbound=s-bound
    pe_predict = (vbg*unbound) + ((np.dot(vec6,j))*v1 *bound) + ((np.dot(vec7,j))*v2 *bound) + \
         ((np.dot(vec8,j))*v3*bound)
    #print(pe_predict.shape)

    diff=pe_predict-p
    return diff
# bound1=([1,1,1,1,1,0,0,0,0.01],[10**15,10**15,10**15,10**15,10**15,10**8,10**8,10**8,1])
res,std=curve_fit(substratebound,data,ydata,p0=[10**4,10**4,10**4,10**4,10**4,100,100,100])
std = np.sqrt(np.diag(std))
print(res)
print(std)
#
# def cost(x):
#     alpha1, alpha2, alpha3, alpha4, alpha5, v1, v2, v3=x
#     vbg=0
#     k=data[0,:]
#     s=data[1,:]
#     p=data[2,:]
#     i=data[3:(len(list1)+3),:]
#     j=data[(len(list1)+3):data.shape[0],:]
#
#     sum = k + s + (np.dot(vec1,i)*alpha1) + (np.dot(vec2,i)*alpha2) + (np.dot(vec3,i)*alpha3) + (np.dot(vec4,i)*alpha4) + (np.dot(vec5,i)*alpha5)
#     bound=0.5 * (sum - ((sum ** 2) - 4 * s * k) ** 0.5)
#
#     unbound=s-bound
#     pe_predict = (vbg*unbound) + ((np.dot(vec6,j))*v1 *bound) + ((np.dot(vec7,j))*v2 *bound) + \
#          ((np.dot(vec8,j))*v3*bound)
#
#     diff=pe_predict-p
#     diff=np.power(diff,2)
#     diff=np.sum(diff)
#     return diff
# res=minimize(cost,[10**4,10**4,10**4,10**4,10**4,10,10,10])
# res=basinhopping(cost,[10**4,10**4,10**4,10**4,10**4,10,10,10])
# print(res)
#
