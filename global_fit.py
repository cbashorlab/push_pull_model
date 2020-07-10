import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import model_alphafit
from scipy.optimize import curve_fit


def onehotencoding(list):
    name = list[0]
    data = pd.read_csv(name)
    factor = pd.Series([name] * len(data['Kinase']), name='Factor')
    data = pd.concat([data, factor], axis=1)
    del list[0]
    for name in list:
        data1 = pd.read_csv(name)
        factor = pd.Series([name] * len(data1['Kinase']), name='Factor')
        data1 = pd.concat([data1, factor], axis=1)
        data = pd.concat([data, data1])
    factor = pd.get_dummies(data.Factor)
    data = data.drop(columns=['Factor'])
    data = pd.concat([data, factor], axis=1)

    vec = pd.get_dummies(['a8v2.csv', 'a8v4.csv', 'a8v9.csv'])
    vec = vec.values

    return data, vec

FileNameList=['a8v2.csv','a8v4.csv','a8v9.csv']
data,vec=onehotencoding(FileNameList)

vec1=vec[0,:]
vec2=vec[1,:]
vec3=vec[2,:]


def substratebound(data,alpha,v1,v2,v3):
    k=data[0,:]
    s=data[1,:]
    i=data[2:5,:]
    sum = k + s + alpha
    pe = (np.dot(vec1,i))*v1 *0.5 * (sum - ((sum ** 2) - 4 * s * k) ** 0.5)+ ((np.dot(vec2,i))*v2 *0.5 * (sum - ((sum ** 2) - 4 * s * k) ** 0.5)) + \
         ((np.dot(vec3,i))*v3*0.5 * (sum - ((sum ** 2) - 4 * s * k) ** 0.5))
  
    return pe


kinase=data['Kinase']
substrate=data['Substrate']
pe=data['PE']
factor=data.iloc[:,4:(5+len(FileNameList))]
kinase=kinase.values
substrate=substrate.values
pe=pe.values
factor=factor.values

xdata=np.vstack([kinase,substrate,factor.T])

popt,pcov=curve_fit(substratebound,xdata,pe,p0=[10**8,10,10,10],bounds=([10,10,10,10],[10**10,10**10,10**10,10**10]))
std = np.sqrt(np.diag(pcov))
print(popt)
print(std)
plt.bar(['Alpha','V1','V2','V3'],popt,yerr=std)
plt.scatter(['Alpha'],10**8,s=100)
plt.scatter(['V1'],10**2,s=100)
plt.scatter(['V2'],10**4,s=100)
plt.scatter(['V3'],10**9,s=100)
plt.yscale('log')
plt.ylabel('Value')
