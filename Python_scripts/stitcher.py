###################################
# THIS SCRIPT STITCHES TOGETHER OUTPUT FROM GHOST INTO A NUMPY ARRAY, SAVES AND PRODUCES PLOTS
# Author: Adrian van Kan
################################### 
import numpy as np
import scipy.io
from matplotlib import pyplot as plt
#import cmocean
import matplotlib as mpl

###############################################################################
######### pyplot settings settings ############################################
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text', usetex=True)

#set font sizes
SMALL_SIZE = 25
MEDIUM_SIZE = 25

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=15)    # legend fontsize
#########################

procs= 3
nx   = 64
ny   = 64
fn   = './ps'

fi_sta   = 6
fi_end   = 7
pre = ''


#########################################################################################################
# READ DATA #
#########################################################################################################
if True:
    for fi in range(fi_sta,fi_end):
       print(fi)
       data=np.zeros((ny,nx))

       ni0=0;
       for i in range(procs):
           f=scipy.io.FortranFile(fn+'.'+str(i).zfill(3)+'.'+pre+str(fi).zfill(3-len(pre))+'.out')
           data_i=f.read_reals(float).reshape((-1,nx,))
           ni=np.size(data_i,0)
           print(np.shape(data_i))
           data[ni0:ni0+ni,:]=data_i
           ni0=ni0+ni

       # OUTPUT fields as *.npy
       np.save(fn+'.stitch.'+pre+str(fi).zfill(3-len(pre))+'.npy',data)
       # OUTPUT fields as *.mat 
#    scipy.io.savemat(fn+'.stitch.'+pre+str(fi).zfill(3-len(pre))+'.mat', {'data': data})
#######################################################################################################

######################################################################################################
# PLOT DATA
######################################################################################################
fn   = 'ps'
Qx   = 1.
Qy   = 1.

for fi in range(fi_sta,fi_end):
   print(fi)
   plt.clf()
   field = np.load(fn+'.stitch.'+pre+str(fi).zfill(3-len(pre))+'.npy')
   plt.pcolormesh(field,cmap='seismic'); plt.show()
