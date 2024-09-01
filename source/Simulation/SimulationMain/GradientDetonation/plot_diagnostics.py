import sys
import yt
import matplotlib.pyplot as plt
ds0 = yt.load('graddet_hdf5_chk_0000')
dd0 = ds0.all_data()
ds1 = yt.load('graddet_hdf5_chk_0001')
dd1 = ds1.all_data()
ds2 = yt.load('graddet_hdf5_chk_0002')
dd2 = ds2.all_data()
ds4 = yt.load('graddet_hdf5_chk_0004')
dd4 = ds4.all_data()


f, ((axtop1,axmid1) ,( axlowmid,axbot)) = plt.subplots(2, 2, sharex=True, figsize=(12,8) )

axtop2 = axtop1.twinx()
axmid2 = axmid1.twinx()

axtop1.plot(dd0['index','x'],dd0['flash','dens'], color='r', linestyle=':')
axtop1.plot(dd1['index','x'],dd1['flash','dens'], color='r', linestyle='-.')
axtop1.plot(dd2['index','x'],dd2['flash','dens'], color='r', linestyle='--')
axtop1.plot(dd4['index','x'],dd4['flash','dens'], color='r', linestyle='-', label=r'$\rho$')
axtop1.set_ylabel('Density (g/cc)', color='red')
axtop1.legend(loc='upper right')
axtop2.plot(dd0['index','x'],dd0['flash','temp'],color='b', linestyle=':')
axtop2.plot(dd1['index','x'],dd1['flash','temp'],color='b', linestyle='-.')
axtop2.plot(dd2['index','x'],dd2['flash','temp'],color='b', linestyle='--')
axtop2.plot(dd4['index','x'],dd4['flash','temp'],color='b', linestyle='-', label='T')
axtop2.set_ylabel('Temperature (K)', color='blue')
axtop2.legend(loc='lower right')


axmid1.plot(dd0['index','x'],dd0['flash','pres'], color='r', linestyle=':')
axmid1.plot(dd1['index','x'],dd1['flash','pres'], color='r', linestyle='-.')
axmid1.plot(dd2['index','x'],dd2['flash','pres'], color='r', linestyle='--')
axmid1.plot(dd4['index','x'],dd4['flash','pres'], color='r', linestyle='-', label='P')
axmid1.set_ylabel('Pressure (erg/cc)', color='red')
axmid2.plot(dd0['index','x'],dd0['flash','velx'],color='b', linestyle=':')
axmid2.plot(dd1['index','x'],dd1['flash','velx'],color='b', linestyle='-.')
axmid2.plot(dd2['index','x'],dd2['flash','velx'],color='b', linestyle='--')
axmid2.plot(dd4['index','x'],dd4['flash','velx'],color='b', linestyle='--', label='vel')
axmid2.set_ylabel('Velocity (cm/s)', color='blue')


axlowmid.plot(dd0['index','x'],dd0['flash','c12 '], color='blue', linestyle=':')
axlowmid.plot(dd0['index','x'],dd0['flash','ne20'], color='green', linestyle=':')
axlowmid.plot(dd0['index','x'],dd0['flash','mg24'], color='red', linestyle=':')
axlowmid.plot(dd0['index','x'],dd0['flash','si28'], color='cyan', linestyle=':')
axlowmid.plot(dd0['index','x'],dd0['flash','ni56'], color='magenta', linestyle=':')

axlowmid.plot(dd1['index','x'],dd1['flash','c12 '], color='blue', linestyle='-.')
axlowmid.plot(dd1['index','x'],dd1['flash','ne20'], color='green', linestyle='-.')
axlowmid.plot(dd1['index','x'],dd1['flash','mg24'], color='red', linestyle='-.')
axlowmid.plot(dd1['index','x'],dd1['flash','si28'], color='cyan', linestyle='-.')
axlowmid.plot(dd1['index','x'],dd1['flash','ni56'], color='magenta', linestyle='-.')

axlowmid.plot(dd2['index','x'],dd2['flash','c12 '], color='blue', linestyle='--')
axlowmid.plot(dd2['index','x'],dd2['flash','ne20'], color='green', linestyle='--')
axlowmid.plot(dd2['index','x'],dd2['flash','mg24'], color='red', linestyle='--')
axlowmid.plot(dd2['index','x'],dd2['flash','si28'], color='cyan', linestyle='--')
axlowmid.plot(dd2['index','x'],dd2['flash','ni56'], color='magenta', linestyle='--')

axlowmid.plot(dd4['index','x'],dd4['flash','c12 '], label='$^{12}$C', color='blue')
axlowmid.plot(dd4['index','x'],dd4['flash','ne20'], label='$^{20}$Ne', color='green')
axlowmid.plot(dd4['index','x'],dd4['flash','mg24'], label='$^{24}$Mg', color='red')
axlowmid.plot(dd4['index','x'],dd4['flash','si28'], label='$^{28}$Si', color='cyan')
axlowmid.plot(dd4['index','x'],dd4['flash','ni56'], label='$^{56}$Ni', color='magenta')
axlowmid.set_ylabel('Fraction')
axlowmid.set_yscale('log')
axlowmid.set_ylim([1e-6,1.5])
axlowmid.legend(loc='best')

x0=dd0['index','x']
x1=dd1['index','x']
x2=dd2['index','x']
x4=dd4['index','x']

dx0 = range(len(x0))
dx1 = range(len(x1))
dx2 = range(len(x2))
dx4 = range(len(x4))
dx0[0] = x0[1]-x0[0]
dx1[0] = x1[1]-x1[0]
dx2[0] = x2[1]-x2[0]
dx4[0] = x4[1]-x4[0]
for i in range(1,len(x0)) :
	dx0[i] = x0[i]-x0[i-1]
for i in range(1,len(x1)) :
	dx1[i] = x1[i]-x1[i-1]
for i in range(1,len(x2)) :
	dx2[i] = x2[i]-x2[i-1]
for i in range(1,len(x4)) :
	dx4[i] = x4[i]-x4[i-1]

axbot.plot(x0,dx0, linestyle=':')
axbot.plot(x1,dx1, linestyle='-.')
axbot.plot(x2,dx2, linestyle='--')
axbot.plot(x4,dx4, linestyle='-')
axbot.set_yscale('log')
axbot.set_ylabel('Cell size (cm)')
axbot.set_xlabel('position (cm)')

#axbot.set_xlim([0,3.6e10])
axbot.set_xlim([0,6e7])

#ax1.set_ylim([1e15,1e25])



plt.show()
