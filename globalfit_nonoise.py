import pandas as pd
import numpy as np
import NoiseModel
from scipy.optimize import curve_fit


def onehotencoding(list,list2):
    vec = pd.get_dummies(list)
    vec = vec.values
    list2.append(list[0])
    vec2 = pd.get_dummies(list2)
    vec2 = vec2.values
    del list2[len(list2) - 1]

    default = list[0]
    data = pd.read_csv(default)
    data = data[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]
    factor_zipper = pd.Series([default] * len(data), name='Zipper')
    factor_velocity = pd.Series(['1'] * len(data), name='Velocity')
    data = pd.concat([data, factor_zipper, factor_velocity], axis=1)

    del list[0]
    for name in list:
        data1 = pd.read_csv(name)
        data1 = data1[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]
        factor_zipper = pd.Series([name] * len(data1), name='Zipper')
        factor_velocity = pd.Series(['1'] * len(data1), name='Velocity')
        data1 = pd.concat([data1, factor_zipper,factor_velocity], axis=1)
        data = pd.concat([data, data1])

    name = list2[0]
    data1 = pd.read_csv(name)
    data1 = data1[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]

    factor_zipper = pd.Series([default] * len(data1), name='Zipper')
    factor_velocity = pd.Series(['2'] * len(data1), name='Velocity')
    data1 = pd.concat([data1, factor_zipper, factor_velocity], axis=1)
    data = pd.concat([data, data1])

    del list2[0]
    velo_index=2
    for name in list2:
        velo_index=velo_index+1
        data1 = pd.read_csv(name)
        data1 = data1[['Flag : Alexa Fluor 680 - Area', 'Myc : Alexa Fluor 750 - Area', 'hospho : PE - Area']]

        factor_zipper = pd.Series([default] * len(data1), name='Zipper')
        factor_velocity = pd.Series([str(velo_index)] * len(data1), name='Velocity')
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
data,vec,vec_va=onehotencoding(list1,list2)

list1=['Globalfit_data/170_127.csv','Globalfit_data/170EEV.csv','Globalfit_data/170EES.csv','Globalfit_data/170EEK.csv','Globalfit_data/170RR.csv']
list2=['Globalfit_data/170L1.csv','Globalfit_data/170L2.csv']
#one hot for zipper
vec1=vec[0,:]
vec2=vec[1,:]
vec3=vec[2,:]
vec4=vec[3,:]
vec5=vec[4,:]
#one hot for va
vec6=vec_va[0,:]
vec7=vec_va[1,:]
vec8=vec_va[2,:]


def substratebound(data,alpha1,alpha2,alpha3,alpha4,alpha5,v1,v2,v3):
    vbg=0
    k=data[0,:]
    s=data[1,:]
    p=data[2,:]
    i=data[3:(len(list1)+3),:]
    j=data[(len(list1)+3):data.shape[0],:]
    sum = k + s + (np.dot(vec1,i)*alpha1) + (np.dot(vec2,i)*alpha2) + (np.dot(vec3,i)*alpha3) + (np.dot(vec4,i)*alpha4) + (np.dot(vec5,i)*alpha5)
    bound=0.5 * (sum - ((sum ** 2) - 4 * s * k) ** 0.5)
    unbound=s-bound
    pe = (vbg*unbound) + ((np.dot(vec6,j))*v1 *bound) + ((np.dot(vec7,j))*v2 *bound) + \
         ((np.dot(vec8,j))*v3*bound)
    pe=pe-p
    return pe



kinase=data['Flag : Alexa Fluor 680 - Area']
substrate=data['Myc : Alexa Fluor 750 - Area']
pe=data['hospho : PE - Area']
factor=data.iloc[:,3:data.shape[1]]
kinase=kinase.values
substrate=substrate.values
pe=pe.values
factor=factor.values


combinedata=np.vstack([kinase,substrate,pe,factor.T])
#Format: Kinase, Substrate, PE, one hot factor
#ydata=combinedata[2,:]
#xdata=np.delete(combinedata,2,0)
ydata=np.zeros(len(combinedata.T))
bound1=([1,1,1,1,1,1,1,1],[10**15,10**15,10**15,10**15,10**15,10**8,10**8,10**8])
popt,pcov=curve_fit(substratebound,combinedata,ydata,p0=[1000,1000,1000,1000,1000,10,10,10])
std = np.sqrt(np.diag(pcov))
print(popt)
print(std)

