import pandas as pd
import numpy as np
import os
#import matplotlib.pyplot as plt
os.chdir('C:\\Users\\Giorgio\\Desktop\\BATCH12_15Gen2020_O3treatement_XR\\')
directory = "C:\\Users\\Giorgio\\Desktop\\BATCH12_15Gen2020_O3treatement_XR\\S11_O30s\\0.4mm_s\\ofet16\\"
file_name = "07_S11_OFET16_VG-5_VD-20_60sON_60sOFF_40kV_100uA.dat"
file = pd.read_csv(directory+file_name,sep='\t', skiprows=11, header=None)
dizi=dict()

##%% #read polynomial data file
with open('polynomial_fit_data_and_mean_current.txt','r') as datafile:
    dizi=dict(eval(datafile.read()))
##%% #delete a row
#dizi.pop(file_name+'_fit_param')   
#dizi.pop(file_name+'_mean_curr')   
##%% 

#polynomial fit of minima
num_cycles = 5
time_cycle = 60 #sec
minima = [file[2][4*(i+1)*time_cycle-1] for i in range(num_cycles)]
times = [file[0][4*(i+1)*time_cycle-1] for i in range(num_cycles)]
fit = list(np.polyfit(times[1:num_cycles],minima[1:num_cycles],2))
fit_short= [float("{:.8e}".format(fit[0])),
            float("{:.8e}".format(fit[1])),
            float("{:.8e}".format(fit[2]))]

#mean current calculation
time_data=file[0]
norm_data=file[2]-(fit[0]*time_data*time_data+fit[1]*time_data+fit[2])

intensities=[float("{:.6}".format(min(norm_data[i*4*time_cycle:(i+1)*4*time_cycle]))) for i in range(2,5)]
mean=        float("{:.6e}".format(np.mean(intensities)))
error=       float("{:.6e}".format((max(intensities)-min(intensities))/2))
mae=         float("{:.6e}".format(np.mean([abs(mean-value) for value in intensities])))

dizi[file_name+'_fit_param']=str(fit_short)+'                 '
dizi[file_name+'_mean_curr']=str([intensities,mean,mae])

##%% #update polynomial data file
with open('polynomial_fit_data_and_mean_current.txt','w') as datafile:
    datafile.write(str(dizi))