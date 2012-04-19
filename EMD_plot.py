from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import EMD

# Signal creation
#base = np.linspace(0,5000,1000)
#a = 50. * np.sin(base/180.) 
#b = 100. * np.sin(base/300.)
#period = (base[1]-base[0])/60.
#int_data = area_data = a + b + 50.*np.random.random(len(base))

#    # Real data
#int_data = np.load('/home/nabobalis/Dropbox/Int_red130.npy')
#area_data = np.load('/home/nabobalis/Dropbox/Area_red130.npy')
int_data = np.load('/home/stuart/Documents/Ellerman_data/IBIS/Int_red130.npy')
area_data = np.load('/home/stuart/Documents/Ellerman_data/IBIS/Area_red130.npy')    
ind = 3801
fin = np.isfinite(int_data[ind,:])
int_data = int_data[ind,:][fin]    
area_data = area_data[ind,:][fin]
time = (np.arange(0,len(int_data))*26.9)/60.
period = 26.9 / 60.

# Calls our EMD!    
imfs_area = EMD.emd(area_data,extrapolation='mirror')
imfs_int = EMD.emd(int_data,extrapolation='mirror')
ncomp = imfs_int.shape[1]
ncomp1 = imfs_area.shape[1]

#if ncomp != ncomp1:
#    print 'End of the line',ncomp,ncomp1
#    raise Exception('Need a better plotting routine')
    
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
    ax1.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
    pltA, = ax1.plot(imfs_area[:,i],'o-',color='r',label="Area")
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
    pltI, = ax2.plot(imfs_int[:,i],'o-',color='b',label="Intensity")        
    j = i + 1
    plt.title("IMF %i"%j)
    if (i+1) == ncomp:
        plt.title("Residual")
    ax1.set_ylabel("Area (pixels)", color='r')
    ax2.set_ylabel("Intensity (A.U.)", color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('r')
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
    ax1.set_xlabel("Frame Number")
#plt.show()

# Signal Recreation
print "-------------------------------------"    
print 'Difference between the signal and the EMD'
print 'Area',np.mean(area_data - np.sum(imfs_area,axis=1))
print 'Intensity',np.mean(int_data - np.sum(imfs_int,axis=1))
print "-------------------------------------"

plt.figure(dpi=80, facecolor='w', edgecolor='k')
ax1 = plt.subplot(111)
ax1.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
ax1.plot(np.sum(imfs_area,axis=1),'ro',label="Area")
ax1.plot(area_data,'r-')   
ax2 = ax1.twinx()
ax2.yaxis.set_major_locator(tick.MaxNLocator(7,symmetric=True))
ax2.plot(np.sum(imfs_int,axis=1),'bo',label="Intensity")
ax2.plot(int_data,'b-')        
ax1.set_ylabel("Area (pixels)", color='r')
ax2.set_ylabel("Intensity (A.U.)", color='b')
ax1.set_xlabel("Frame Number") 
for tl in ax1.get_yticklabels():
    tl.set_color('r')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
#plt.show()

# FFT Analysis
plt.figure(figsize=(16, 10), dpi=80, facecolor='w', edgecolor='k')
plt.subplots_adjust(wspace = 0.4)
a = np.ceil(ncomp/2.)
n= len(imfs_area)
freq = np.arange(n/2)/(n*period)

for i in np.arange(0,ncomp-1):
    ax1 = plt.subplot(2,a,i+1)
    fft_a = np.abs(np.fft.fft(imfs_area[:,i])[0:n/2])#**2
    fft_i = np.abs(np.fft.fft(imfs_int[:,i])[0:n/2])#**2        
    pltA, = ax1.plot((1./freq),fft_a,'o-',color='r',label="Area")
    ax2 = ax1.twinx()
    pltI, = ax2.plot((1./freq),fft_i,'o-',color='b',label="Intensity")        
    j = i + 1
    plt.title("IMF %i"%j)
    ax1.set_xlabel("Period [Minutes]")
    ax1.set_ylabel("Power [A.U.]", color='r')
    ax2.set_ylabel("Power [A.U.]", color='b')
    ax1.axis(xmin=0,xmax=15)
    for tl in ax1.get_yticklabels():
        tl.set_color('r')
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
plt.show()