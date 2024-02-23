import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import scipy
from scipy.io import FortranFile

def interpolate(field,res_old,Qs_old,res_new,Qs_new):
   nx_old, ny_old = res_old
   Qx_old, Qy_old = Qs_old
   x = np.linspace(0, 2*np.pi/Qx_old, nx_old)        # the grid is an outer product
   y = np.linspace(0, 2*np.pi/Qy_old, ny_old)        # of x and y arrays

   func = RectBivariateSpline(x, y, field, s=0)
   nx, ny = res_new
   Qx, Qy = Qs_new
   xnew = np.linspace(0,2*np.pi/Qx,nx)
   ynew = np.linspace(0,2*np.pi/Qy,ny)
   field_new = func(xnew, ynew)
   return field_new

def read_field_from_file(fn,fi,procs,nx,ny):
     data=np.zeros((ny,nx))
     ni0=0; pre = ''
     for i in range(procs):
       f=scipy.io.FortranFile(fn+'.'+str(i).zfill(3)+'.'+pre+str(fi).zfill(3-len(pre))+'.out')
       data_i=f.read_reals(float).reshape((-1,nx))
       ni=np.size(data_i,0)
       data[ni0:ni0+ni,:]=data_i
       ni0=ni0+ni
     return data

#THE PROGRAMS STARTS
field = read_field_from_file(fn='./ps',fi=5,procs=4,nx=32,ny=32)
Qs_old = [1,1]
res_old = [np.shape(field)[0],np.shape(field)[1]]
nx,ny = res_old
#x = np.linspace(0,2*np.pi/Qs_old[0],nx); y = np.linspace(0,2*np.pi/Qs_old[1],ny)
#X,Y = np.meshgrid(x,y)

## CHANGE RESOLUTION HERE
fact_x = 2
fact_y = 2
Qs_new  = [1,1]
res_new = [fact_x*res_old[0],fact_y*res_old[1]]
field_scaled_up = interpolate(field,res_old,Qs_old,res_new,Qs_new)
if False: #True for graphical check
  plt.figure(0)
  plt.pcolormesh(field)
  plt.figure(1)
  plt.pcolormesh(field_scaled_up)
  plt.show()

## NOW SAVE RESULT TO FILE with arbitrary number of CPUs according to GHOST's parallelization
nprocs_new = 3  #NEW NUMBER OF CPUS
fi_out     = 6 #index of regridded field output
n1 = 0;  n2 = res_new[1]-1;
print('n1,n2=',n1,n2)
iwork1 = (n2-n1+1)//nprocs_new
iwork2 = np.mod(n2-n1+1,nprocs_new)
print('iwork1,iwork2=',iwork1,iwork2)
for irank in range(nprocs_new):
  jsta = int(irank*iwork1+n1+min(irank,iwork2))
  jend = int(jsta+iwork1-1)
  if (iwork2 > irank): 
    jend = jend+1
    print('nonono')
  print('jsta,jend=',jsta,jend)
  f = FortranFile('ps.'+str(irank).zfill(3)+'.'+str(fi_out).zfill(3)+'.out', 'w')
  f.write_record(field_scaled_up[jsta:jend+1,:])
  f.close()


