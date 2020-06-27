import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
import model_alphafit
#import random

matrix_noise = pd.read_csv('HANoise.csv')

matrix_noise = matrix_noise[matrix_noise['CD8 : APC - Area'] > 0]
matrix_noise = matrix_noise[matrix_noise['GFP - Area'] > 0]
myc = matrix_noise['CD8 : APC - Area']
GFP = matrix_noise['GFP - Area']

list=model_alphafit.noiselevel(myc, GFP, 100, 100)



def alphaglobalfit(File_name,kinase_name,substrate_name,PE_name,va,noise_k,noise_s,noise_pe):
    data=pd.read_csv(File_name)
    data = data[data[kinase_name] > 0]
    data = data[data[substrate_name] > 0]
    data = data[data[PE_name] > 0]
    # positions = []
    # for index in range(1):
    #     positions.append(random.randint(1, 6000))
    # data = data.iloc[positions,:]
    kinase=data[kinase_name]
    substrate=data[substrate_name]
    PE=data[PE_name]
    mean_alpha, prob_alpha = model_alphafit.multicell(noise_k,noise_s,noise_pe,kinase,substrate,PE,va,100)
    return mean_alpha,prob_alpha



va=np.arange(1,100,1)
#mean_alphaV, prob_alphaV=alphaglobalfit('170EEV.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
#mean_alphaL, prob_alphaL=alphaglobalfit('170EEL.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
#mean_alphaN, prob_alphaN=alphaglobalfit('170EEN.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
#mean_alphaS, prob_alphaS=alphaglobalfit('170EES.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
#mean_alphaK, prob_alphaK=alphaglobalfit('170EEK.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
combine=np.array([0,0,0])
for value in va:
    mean_alphaR, prob_alphaR=alphaglobalfit('170EEV.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',value,list,list,list)
    va_store=np.ones(len(prob_alphaR))*value
    total=np.vstack([mean_alphaR,va_store,prob_alphaR])
    total=total.T
    combine=np.vstack([combine,total])
np.savetxt('Result.csv',combine)
# plt.plot(mean_alphaR,prob_alphaR)
# plt.xscale('log')

# mean_alphaL, prob_alphaL=vaglobalfit('168.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
# mean_alphaN, prob_alphaN=vaglobalfit('169.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
# mean_alphaS, prob_alphaS=vaglobalfit('170.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
# mean_alphaK, prob_alphaK=vaglobalfit('171.csv','Flag : Alexa Fluor 680 - Area','Myc : Alexa Fluor 750 - Area','hospho : PE - Area',va,list,list,list)
#

# plt.figure(1)
# plt.plot(mean_alphaV,prob_alphaV)
# plt.plot(mean_alphaL,prob_alphaL)
# plt.plot(mean_alphaN,prob_alphaN)
# plt.plot(mean_alphaS,prob_alphaS)
# plt.plot(mean_alphaK,prob_alphaK)
# plt.plot(mean_alphaR,prob_alphaR)
# plt.xlabel('Alpha')
# plt.ylabel('Normalized Probablity')
# plt.title('Alpha Probablity')
# plt.xscale('log')

# max_index=np.argmax(prob_alphaV)
# va_max167=mean_alphaV[max_index]
# max_index=np.argmax(prob_alphaL)
# va_max168=mean_alphaL[max_index]
# max_index=np.argmax(prob_alphaN)
# va_max169=mean_alphaN[max_index]
# max_index=np.argmax(prob_alphaS)
# va_max170=mean_alphaS[max_index]
# max_index=np.argmax(prob_alphaK)
# va_max171=mean_alphaK[max_index]
# max_index=np.argmax(prob_alphaR)
# va_max172=mean_alphaR[max_index]
#
# Name=['EEV','EEL','EEN','EES','EEK','RR']
#
# height=np.array([va_max167,va_max168,va_max169,va_max170,va_max171,va_max172])
# plt.figure(2)
# position=np.array([1,2,3,4,5,6])
# plt.bar(position,height)
# plt.xticks(position,Name)
# plt.ylabel('Alpha')
# toc=time.perf_counter()
# print(toc-tic)