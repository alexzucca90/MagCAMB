#Loading libraries.. Maybe you have to download them..
import numpy as np
import matplotlib as mpl
import pylab as pil
import matplotlib.pyplot as pl
# this you don't need now.. It's for multiple plots..
#from matplotlib import gridspec


#Thicker axes
MP_LINEWIDTH = 2
pil.rc('axes', linewidth = MP_LINEWIDTH)

label_size = 10
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

# loading data
# change the names between ' '
#data1 = np.loadtxt('scalar_spectrum.txt')
#data2 = np.loadtxt('tensor_spectrum.txt')
#data3 = np.loadtxt('Bmodes_tensCls.dat')
#data4 = np.loadtxt('Bmodes_lensedCls.dat')
data1 = np.loadtxt('SPT_combined.dat')
#data2 = np.loadtxt('B_modes.dat')
data3 = np.loadtxt('VecComp-29_vecCls.dat')
data4 = np.loadtxt('VecComp-25_vecCls.dat')
data5 = np.loadtxt('test_lensedCls.dat')
data6 = np.loadtxt('VecComp-20_vecCls.dat')
data7 = np.loadtxt('VecComp-15_vecCls.dat')

#pl.plot for linear scale, pl.loglog for log log scale, pl.semilogyfor log in the y axis only
#pl.loglog(data1[:,0], data1[:,0]*(data1[:,0]+1)*data1[:,1]/6.28, 'red', label = r"scalar")
#pl.loglog(data2[:,0], data2[:,0]*(data2[:,0]+1)*data2[:,1]/6.28, 'blue', label = r"tensor")
#pl.loglog(data3[:,0], data3[:,3], 'green', label = r"lensed")
#pl.loglog(data4[:,0], data4[:,3], 'black', label = r"tensor")
#pl.semilogx(data2[:,0], data2[:,1], 'black', label = r" clean ")
#pl.semilogx(data2[:,0], data2[:,2], 'red', label = r" + dust & noise")
pl.semilogx(data1[:,0], data1[:,1], 'Black', marker= 's', linestyle = "", label = "SPT combined")
pl.errorbar(data1[:,0], data1[:,1], yerr=data1[:,2],linestyle = "", color = 'Black')
#pl.semilogx(data1[:,0], data1[:,2], 'green', marker= 'p', linestyle = "", label = "SPT 95x150")
#pl.semilogx(data1[:,0], data1[:,3], 'yellow', marker= 'h', linestyle = "", label = "SPT 150x150")
#pl.semilogx(data3[:,0], data3[:,1], 'Black', linestyle = "--", label = r"Tens+Lens primary")
pl.semilogx(data5[:,0], data5[:,3]*6.28/(data5[:,0]+1)*1000, 'Grey', linestyle = "-", linewidth="2", label = r"Lensing B-modes")
pl.semilogx(data3[:,0], data3[:,3]*6.28/(data3[:,0]+1)*1000, 'Red', linestyle = "-", linewidth="2", label = r"$n_B=-2.9$")
pl.semilogx(data4[:,0], data4[:,3]*6.28/(data4[:,0]+1)*1000, 'Blue', linestyle = "-", linewidth ="2", label = r"$n_B=-2.5$")
pl.semilogx(data6[:,0], data6[:,3]*6.28/(data6[:,0]+1)*1000, 'Green', linestyle = "-", linewidth ="2", label = r"$n_B=-2.0$")
pl.semilogx(data7[:,0], data7[:,3]*6.28/(data7[:,0]+1)*1000, 'Orange', linestyle = "-", linewidth ="2", label = r"$n_B=-1.5$")
#pl.semilogx(data4[:,0], data4[:,3]*6.28/(data4[:,0]+1)*1000, 'Red', linestyle = "--", label = r"Mag Vec Comp B-modes")
#pl.semilogx(data6[:,0], data6[:,3]*6.28/(data6[:,0]+1)*1000, 'Red', linestyle = "-", label = r"Mag Ten Pass B-modes")
#pl.semilogx(data5[:,0], (data4[:,3]+data5[:,3])*1000*6.28/(data5[:,0]+1), 'Red', linestyle = "-", label = r"SUM-Bmodes")
#pl.semilogx(data5[:,0], data5[:,1]*6.28/(data5[:,0]+1), 'Red', linestyle = "-.", label = r"Tens+Lens primary")
# Labels
pil.xlabel(r'$l$', fontsize = 25)
pil.ylabel(r'$l C_l^{BB} \, [10^{-3} (\mu K^2)]$', fontsize = 25)
pl.xlim(50,5000)
pl.ylim(-1,2.0)
# where you want to put the legend
pl.legend(loc='upper left', shadow=True)

# Thicker Markers on the axes
pil.tick_params(axis='both',color='black',length=5,width=2)

#Lowering the x-axes labels
pil.tick_params(axis='x', pad=8)
pil.tick_params(axis='y', pad=8)
pil.tick_params(axis='both', which='major', labelsize=15)
# Saving
pil.savefig('SPT_MagInd.pdf', bbox_inches='tight')
