import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
import matplotlib.pyplot as plt

diz_directory= dict()
diz_directory[0]="C:\\Users\\Giorgio\\Desktop\\BATCH12_15Gen2020_O3treatement_XR\\"

file_names= dict()
################################ S11 
file_names[11_0.4_13_0]='02_S11_OFET13_trfSat_preXR'   #35..25..80..89
file_names[11_0.4_13_1]='02_S11_OFET13_trfSat_postXR'  #35..25..80..89
file_names[11_0.4_16_0]='02_S11_OFET16_trfSat_preXR'   #35..25..85..95
file_names[11_0.4_16_1]='02_S11_OFET16_trfSat_postXR'  #35..25..85..95


directory = diz_directory[0]
os.chdir(directory)

key=11_0.4_16_1

file_name = file_names[key]
extension = '.DAT'

file = pd.DataFrame
file = pd.read_csv(diz_directory[0]+file_name+extension,sep='\t', skiprows=11, header=None)
file.columns = ['Time (ms)','Vdrain_(V)','Idrain_(A)','Vgate_(V)','Igate_(A)'
                ]

#### Constants
e0 = 8.854187e-12 *1e-2  # F/cm
er = 3 # Parylene
W = 0
L = 40 *1e-4 # cm
d = 260 *1e-7 #cm  Parylene 
# s11 260 s12 252 s13 244 s14 240 s15 236
Vd = -20 #V

#if '0.4'==str(key)[2:5]:
#    W = 35 *1e-1 # cm 0.4mm/s
#elif '0.8'==str(key)[2:5]:
#    W = 42 *1e-1 # cm 0.8mm/s
W = 35 *1e-1                                                  ### verify
    
file['Idrain_abs_(A)'] = file['Idrain_(A)'].abs()
file['sqrt_Idrain_(A^1/2)'] = np.sqrt((file['Idrain_abs_(A)']))
file['log_Idrain'] = np.log10((file['Idrain_abs_(A)'])) 


######## ON/OFF ratio
len_interval=5                                  ####### DEFAULT =10
len_col =len(file['Idrain_abs_(A)'])

#off_array = np.array(file['Idrain_abs_(A)'][0:len_interval])   
off_array = np.array(file['Idrain_abs_(A)'][(len_col-len_interval):len_col])   
off = np.mean(off_array)

on_array = np.array(file['Idrain_abs_(A)'][int(len_col/2)-1])
on = float(on_array)

on_off_ratio = on/off

############# max mobility

def lin_fit_mob(x,A,B):
    C=(W*e0*er)/(2*L*d)
#    A=abs(A)
    return -np.sqrt(C*A)*x+np.sqrt(C*A)*B

u_lim_mob = file[file['sqrt_Idrain_(A^1/2)']>0.0035].index[0]
b_lim_mob = file[file['sqrt_Idrain_(A^1/2)']>0.0025].index[0]

lf_mob_opt, lf_mob_cov = curve_fit(lin_fit_mob,
                                   file['Vgate_(V)'][b_lim_mob:u_lim_mob],
                                   file['sqrt_Idrain_(A^1/2)'][b_lim_mob:u_lim_mob],
                                   bounds=([0,-np.inf],[np.inf,np.inf]))

file['sqrt_Idrain_lin_fit'] = lin_fit_mob(file['Vgate_(V)'],lf_mob_opt[0],lf_mob_opt[1])

############# subthreshold slope
base_line=np.log10(off)
#file = file.iloc[::-1].reset_index(drop=True)    ####### DEFAULT COMMENTED
def lin_fit_thr(x,A,B):
    return -(1/A)*x+(1/A)*B+base_line

u_lim_thr = file[file['log_Idrain']>-8.5].index[0]   
b_lim_thr = file[file['log_Idrain']>-9.5].index[2]   

lf_thr_opt, lf_thr_cov = curve_fit(lin_fit_thr,
                                   file['Vgate_(V)'][b_lim_thr:u_lim_thr],
                                   file['log_Idrain'][b_lim_thr:u_lim_thr],
                                   bounds=([0,-np.inf],[np.inf,np.inf]))

file['log10_Idrain_lin_fit'] = lin_fit_thr(file['Vgate_(V)'],lf_thr_opt[0],lf_thr_opt[1])

#file = file.iloc[::-1].reset_index(drop=True)    ####### DEFAULT COMMENTED


##plt.plot(file['Vgate_(V)'],file['Idrain_(A)'])
fig,axs = plt.subplots(nrows=2,figsize=[4.8, 7.2])
axs[0].plot(file['Vgate_(V)'],lin_fit_mob(file['Vgate_(V)'],lf_mob_opt[0],lf_mob_opt[1]),'r', linewidth = 1)
axs[1].plot(file['Vgate_(V)'],lin_fit_thr(file['Vgate_(V)'],lf_thr_opt[0],lf_thr_opt[1]),'r', linewidth = 1)
#axs[0].plot([file['Vgate_(V)'][b_lim_mob],file['Vgate_(V)'][u_lim_mob]],[file['sqrt_Idrain_(A^1/2)'][b_lim_mob],file['sqrt_Idrain_(A^1/2)'][u_lim_mob]])
axs[1].plot([file['Vgate_(V)'][b_lim_thr],file['Vgate_(V)'][u_lim_thr]],[file['log_Idrain'][b_lim_thr],file['log_Idrain'][u_lim_thr]])
axs[0].plot(file['Vgate_(V)'],file['sqrt_Idrain_(A^1/2)'], 'k', linewidth = 1)
axs[1].plot(file['Vgate_(V)'],file['log_Idrain'],'k', linewidth = 1)
axs[1].plot(file['Vgate_(V)'][0:len_interval],file['log_Idrain'][0:len_interval],'b')
axs[0].grid()
axs[1].grid()
axs[0].set_ylabel('sqrt_Idrain')
axs[1].set_ylabel('log_Idrain')
#axs[0].set_xlim(-20,15)
#axs[1].set_xlim(-10,30)
#axs[1].set_ylim(-11,-4)


#### prinf datafile.txt
with open(file_name+"_plus"+'.txt','w') as datafile:
    datafile.write(file.to_string(index= False ,   decimal = ','))


#directory = "C:\\Users\\Giorgio\\Desktop\\BATCH11_17gen2020_TOPGATE"
#os.chdir(directory)


param_file = pd.DataFrame()

try:
    with open('transfer_params.txt','r') as transfer_file:
        param_file= pd.read_csv(transfer_file, header=0, delim_whitespace=True)
#        print(param_file)
except IOError:
    with open('transfer_params.txt','w') as transfer_file:
    
        param_file['a']=[file_name]
        param_file['b']=[0]
        param_file['c']=[0]
        param_file['d']=[0]
        param_file['e']=[0]
        param_file['f']=[0]
        param_file['g']=[0]
        param_file['h']=[0]
        param_file['i']=[0]
        param_file['j']=[0]

    
col_labels =  ['file_name',
               'max_mob_(cm^2/Vs)',
               'mob_err',
               'Vt_(V)',
               'Vt_err',
               'subt.sl_(V/dec)',
               'sub_err',
               'Von_(V)',
               'Von_err',
               'ON/OFF_ratio']    

param_file.columns = col_labels

for name in param_file['file_name']:
    if file_name==name:
        num = param_file.index[param_file['file_name']==file_name].item()
        param_file.loc[num:num,'file_name':'ON/OFF_ratio']= \
        (file_name,
         lf_mob_opt[0],
         np.sqrt(lf_mob_cov[0][0]),
         lf_mob_opt[1],
         np.sqrt(lf_mob_cov[1][1]),
         lf_thr_opt[0],
         np.sqrt(lf_thr_cov[0][0]),
         lf_thr_opt[1],
         np.sqrt(lf_thr_cov[1][1]),
         on_off_ratio)
        break
else:
    new = pd.DataFrame(
      [[file_name,
      lf_mob_opt[0],
      np.sqrt(lf_mob_cov[0][0]),
      lf_mob_opt[1],
      np.sqrt(lf_mob_cov[1][1]),
      lf_thr_opt[0],
      np.sqrt(lf_thr_cov[0][0]),
      lf_thr_opt[1],
      np.sqrt(lf_thr_cov[1][1]),
      on_off_ratio]],
      columns=col_labels)
    param_file = pd.concat([param_file,new], ignore_index=True)

pd.set_option('display.max_colwidth',1000)
with open('transfer_params.txt','w') as transfer_file:
    transfer_file.write(param_file.to_string(index= False))
    
print(param_file[['file_name']])
    