import pandas as pd
import numpy as np
from scipy.optimize import minimize



data=pd.read_csv('PTPN3_V2.csv')
data=data[data['Flag : Alexa Fluor 680 - Area']>0]
data=data[data['Myc : Alexa Fluor 750 - Area']>0]
data=data[data['hospho : PE - Area']>0]
data=data[data['CD8 : APC - Area']>0]
kinase=data['Flag : Alexa Fluor 680 - Area']
substrate=data['Myc : Alexa Fluor 750 - Area']
pe=data['hospho : PE - Area']
phosphatase=data['CD8 : APC - Area']

kinase=kinase.values
substrate=substrate.values
pe=pe.values
phosphatase=phosphatase.values


data=np.vstack([kinase,substrate,phosphatase,pe])
#Format: Kinase, Substrate, PE, one hot factor
#ydata=combinedata[2,:]
#xdata=np.delete(combinedata,2,0)
ydata=np.zeros(len(data.T))

def costfunc(input):
    alpha,alpha2,v_k,v_p=input
    vbg=0
    k=data[0,:]
    s=data[1,:]
    p=data[2,:]
    pe=data[3,:]
    kinase=np.ones(k.shape)
    phosphatase=np.ones(k.shape)
    old_substrate=np.ones(k.shape)
    substrate=s/(kinase/alpha+phosphatase/alpha2+1)
    loop=0
    while abs(np.sum(substrate)-np.sum(old_substrate))>0.01:
        loop=loop+1
        old_substrate=substrate
        kinase=k/(1+(substrate/alpha))
        phosphatase=p/(1+(substrate/alpha2))
        substrate = s / (kinase / alpha + phosphatase / alpha2 + 1)
    predict=s*(v_k*kinase/alpha+vbg)/(v_k*kinase/alpha+vbg+1+v_p*phosphatase/alpha2)
    cost=predict-pe
    cost=np.power(cost,2)
    cost=np.sum(cost)
    print(loop)
    return cost
bound1=((0,None),(0,None),(0,None),(0,None))
# res=minimize(costfunc,[1000,1000,10,10],method='L-BFGS-B',bounds=bound1)
# print(res)
res=minimize(costfunc,[1000,1000,10,10],method='Nelder-Mead')
print(res)