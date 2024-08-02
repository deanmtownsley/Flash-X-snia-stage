import sys
import yt
import matplotlib.pyplot as plt
ds1 = yt.load(sys.argv[1])
dd1 = ds1.all_data()


#ax1= plt.subplot2grid((1,1),(0,0))
f, (axtop1,axmid1,axlowmid,axbot) = plt.subplots(4, sharex=True, figsize=(8,12) )

axtop2 = axtop1.twinx()
axmid2 = axmid1.twinx()

axtop1.plot(dd1['index','x'],dd1['flash','dens'])
#axtop1.plot(dd1['index','r'],dd1['flash','dens'], 'o')
#axtop1.set_yscale('log')
axtop1.set_ylabel('Density (g/cc)')
axtop2.plot(dd1['index','x'],dd1['flash','temp'],color='r')
axtop2.set_ylabel('Temperature (K)')


axmid1.plot(dd1['index','x'],dd1['flash','pres'])
axmid1.set_ylabel('Pressure (erg/cc)')
axmid2.plot(dd1['index','x'],dd1['flash','velx'],color='r')
axmid2.set_ylabel('Velocity (cm/s)')

axlowmid.plot(dd1['index','x'],dd1['flash','c12 '], label='c12')
axlowmid.plot(dd1['index','x'],dd1['flash','ne20'], label='ne20')
axlowmid.plot(dd1['index','x'],dd1['flash','mg24'], label='mg24')
axlowmid.plot(dd1['index','x'],dd1['flash','si28'], label='si28')
axlowmid.plot(dd1['index','x'],dd1['flash','ni56'], label='ni56')
axlowmid.set_ylabel('Fraction')
axlowmid.set_yscale('log')
axlowmid.set_ylim([1e-6,1.5])
axlowmid.legend(loc='best')

x=dd1['index','x']

dx = range(len(x))
dx[0] = x[1]-x[0]
for i in range(1,len(x)) :
	dx[i] = x[i]-x[i-1]

axbot.plot(x,dx)
axbot.set_yscale('log')
axbot.set_ylabel('Cell size (cm)')
axbot.set_xlabel('position (cm)')

#axbot.set_xlim([0,3.6e10])
axbot.set_xlim([0,1e8])

#ax1.set_ylim([1e15,1e25])



plt.show()
