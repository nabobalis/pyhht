
from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from pyhht.EMD import emd

# Signal creation
base = np.linspace(0,5000,1000)
a = 50. * np.sin(base/18.)
b = 100. * np.sin(base/30.)
period = (base[1]-base[0])/60.
int_data = area_data = a + b #+ 5.*np.random.random(len(base))
time = np.arange(0,len(base)*period,period)

# Calls our EMD!
imfs_area = emd(area_data)
imfs_int = emd(int_data)
ncomp = imfs_int.shape[1]
ncomp1 = imfs_area.shape[1]

print "-------------------------------------"
print "Area:",imfs_area.shape
print "Intensity:",imfs_int.shape
print "-------------------------------------"

# Sexy image plotting
plt.figure(figsize=(16, 10), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(wspace = 0.4)
a = np.ceil(ncomp/2.)
for i in np.arange(0,ncomp):
    ax1 = plt.subplot(2,a,i+1)
    try:
        ax1.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
        pltA, = ax1.plot(time,imfs_area[:,i],'-',color='r',label="Area")
    except IndexError:
        pass
    ax2 = ax1.twinx()
    try:
        ax2.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
        pltI, = ax2.plot(time,imfs_int[:,i],'-',color='b',label="Intensity")
    except IndexError:
        pass
    j = i + 1
    plt.title("IMF %i"%j)
    if (i+1) == ncomp:
        plt.title("Residual")
    ax1.set_ylabel("Area [pixels]", color='r')
    ax2.set_ylabel("Intensity [A.U.]", color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('r')
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
#    ax1.set_xlabel("Frame Number")
    ax1.set_xlabel("Time [Minutes]")
plt.show()

# Signal Recreation
print "-------------------------------------"
print 'Difference between the signal and the EMD'
print 'Area',np.mean(area_data - np.sum(imfs_area,axis=1))
print 'Intensity',np.mean(int_data - np.sum(imfs_int,axis=1))
print "-------------------------------------"

plt.figure(dpi=80, facecolor='w', edgecolor='k')
ax1 = plt.subplot(111)
ax1.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
ax1.plot(np.sum(imfs_area,axis=1),'o',color='r',label="Area")
ax1.plot(area_data,'r-')
ax2 = ax1.twinx()
ax2.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
ax2.plot(np.sum(imfs_int,axis=1),'o',color='b',label="Intensity")
ax2.plot(int_data,'b-')
ax1.set_ylabel("Area [pixels]", color='r')
ax2.set_ylabel("Intensity [A.U.]", color='b')
#ax1.set_xlabel("Frame Number")
ax1.set_xlabel("Time [Minutes]")
for tl in ax1.get_yticklabels():
    tl.set_color('r')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
plt.show()

# FFT Analysis
plt.figure(figsize=(16, 10), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(wspace = 0.4)
a = np.ceil(ncomp/2.)
n= len(imfs_area)
freq = np.arange(n/2)/(n*period)

for i in np.arange(0,ncomp):
    ax1 = plt.subplot(2,a,i+1)
    fft_a = np.abs(np.fft.fft(imfs_area[:,i])[0:n/2])#**2
    fft_i = np.abs(np.fft.fft(imfs_int[:,i])[0:n/2])#**2
    pltA, = ax1.plot((1./freq),fft_a,'-',color='r',label="Area")
    ax2 = ax1.twinx()
    pltI, = ax2.plot((1./freq),fft_i,'-',color='b',label="Intensity")
    j = i + 1
    plt.title("IMF %i"%j)
    ax1.set_xlabel("Period [Minutes]")
    ax1.set_ylabel("Power [A.U.]", color='r')
    ax2.set_ylabel("Power [A.U.]", color='b')
    ax1.axis(xmin=0,xmax=60)
    for tl in ax1.get_yticklabels():
        tl.set_color('r')
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
plt.show()