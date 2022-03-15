# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 02:20:21 2022

@author: tynes
"""

import numpy as np
import matplotlib.pyplot as plt

from STGM2_ import sem2d_read_fault
from STGM2_ import get_snapshot_data
from STGM2_ import get_asv_data
from STGM2_ import get_max_contour
from STGM2_ import plot_1
from STGM2_ import plot_2
from STGM2_ import plot_2_xlab
from STGM2_ import plot_contour_v
from STGM2_ import plot_max_contour
from STGM2_ import plot_xmaxtmax_of_dBarr
from STGM2_ import plot_xmaxtmax_of_DcBarr
from STGM2_ import one_anim

from STGM2_ import get_path_STGM2
from STGM2_ import get_Par_data
from STGM2_ import convertir_Par_data
from STGM2_ import all_,all__
from STGM2_ import print_,print__
from STGM2_ import rsf_mu





#------------------------------------------------------------------------------------
#                                '''initialisation '''
#------------------------------------------------------------------------------------
#
#
models_path=get_path_STGM2()

#
#
#
#------------------------------------------------------------------------------------
#
#
m=0
#
#

#------------------------------------------------------------------------------------
#                                '''debut programme '''
#------------------------------------------------------------------------------------
#  
#
D=get_Par_data((models_path[m]))
#print__(D)
D=convertir_Par_data(D)
print_(D)
#
#
#
#all_(models_path,ref='contour')
#all_(models_path,ref='max_contour')
#all__(models_path,ref1='contour',ref2='max_contour')



plot_xmaxtmax_of_dBarr(models_path)  #test_dBarr
#plot_xmaxtmax_of_DcBarr(models_path)#test_DcBarr
#
#

# 
#print(rsf_mu(D)['Tt']) 

#--------------------------------------------------------------------------------------
#all_(models_path,ref='Rup_speed')
#all_(models_path,ref='Rup_acc')
#all_(models_path,ref='G0')
#all_(models_path,ref='Peak_sliprate')
#all_(models_path,ref='slip')
#all_(models_path,ref='stress_drop')
#all_(models_path,ref='sliprate')
#all_(models_path,ref='stress_change')
#
#
#all__(models_path,'sliprate','stress_change')
#all__(models_path,'Rup_speed','Rup_acc')
#all__(models_path,'Rup_speed','G0')
#all__(models_path,'Rup_speed','Peak_sliprate')
#all__(models_path,'Rup_speed','slip')
#all__(models_path,'Rup_speed','stress_drop')
#

#data = sem2d_read_fault(models_path[m],fault_name="Flt01")
#x,sliprate,stress_change=get_snapshot_data(models_path[m])
#xx,Rup_t,Heal_t,Rup_speed,Rup_acc,Gc,G0,Peak_sliprate,slip,stress_drop=get_asv_data(models_path[m])
#
#
#plot_1(X=xx,y=stress_drop,ref='slip',titre=models_path[m])
#plot_1(X=x,y=sliprate,ref='sliprate',titre=models_path[m])
#plot_1(X=x,y=stress_drop,ref='sliprate',titre=models_path[m])
#plot_2(X=x,y1=sliprate,y2=stress_change,ref1='sliprate',ref2='stress_change')
#plot_2(X=xx,y1=Gc,y2=G0,ref1='Gc',ref2='G0')
#plot_2(X=xx,y1=Peak_sliprate,y2=stress_drop,ref1='Peak_sliprate',ref2='stress_drop')
#plot_2(X=xx,y1=Heal_t,y2=Rup_t,ref1='Heal',ref2='Rup')
#plot_contour_v(data)
#
#----animations----------------------------------------------------------------------
#
#one_anim(data,fr=400,ref='v')
#one_anim(data,fr=400,ref='st')
#one_anim(data,fr=400,ref='sn')
#one_anim(data,fr=400,ref='mu')
#one_anim(data,fr=400,ref='d')
#
#
#plot_max_contour(data)

#
#
#

# 
# 