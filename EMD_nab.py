"""
To do:
    Make it work
    Tests / Examples (same in some literature)
    Instanous freequncuies
    Hilbert-Huang Transform
"""
from __future__ import division
import numpy as np
from sklearn.svm import SVR
from scipy import interpolate
import matplotlib.pyplot as plt

__all__ = 'emd'

def peak(data):

    #############################################
    # Defines the envelopes
    #############################################
    min_env = np.zeros([len(data)],dtype=bool)
    max_env = min_env.copy()
    #############################################
    # Finds the extrema of the data
    #############################################
    max_env = np.logical_and(
            np.r_[True, data[1:] > data[:-1]],
            np.r_[data[:-1] > data[1:], True])
    min_env = np.logical_and(
            np.r_[True, data[1:] < data[:-1]],
            np.r_[data[:-1] < data[1:], True])
    #############################################
    # Fixes some issues with the above finding routine
    #############################################
#    if max_env[0] > min_env[0]:
#        min_env[0] =  False
#    else:
#        max_env[0] =  False
    if min_env[-1] > max_env[-1]:
        min_env[-1] =  False
    else:
        max_env[-1] =  False
    min_env = min_env.nonzero()[0]
    max_env = max_env.nonzero()[0]
    min_env = min_env[min_env.nonzero()[0]]
    max_env = max_env[max_env.nonzero()[0]]

    if min_env[0] > max_env[0]:
        first = max_env[0]
    else:
        first = min_env[0]

    if min_env[-1] > max_env[-1]:
        last = max_env[-1]
    else:
        last = min_env[-1]

    return min_env,max_env,first,last

def peak_para(data):

    ypoints = []
    xpoints = []
    points = []
    for i in xrange(1,data.size-1):
        y = np.array([data[i-1],data[i],data[i+1]])
        val = np.array([[0.5, -1, 0.5],[-2.5,4,-1.5],[3,-3,1]])
        a,b,c = np.dot(val,y)
        if a == 0:
            pass
        tp = - b / (2*a)
        if tp <= 2.5 and tp >= 1.5:
            xpoints.append(i-2+tp)
#            ypoints.append(-b*(tp)/2 + c) # WRONG BUT FROM THE PAPER
            ypoints.append(a*(2)**2 + b*(2)+c) # "SOLVES" A PARABOLIC EQN FOR Y - SEEMS TO BE WRITE
            if a < 0:
                points.append('max')
            if a > 0:
                points.append('min')
    return xpoints,ypoints,points


def predict(data):
    svr_rbf = SVR(kernel='rbf', C=1e3, gamma=0.1)
    X = np.linspace(0,data.size,data.size)[:,None]
    start = np.linspace(-data.size*0.1,0,100)[:,None]
    end = np.linspace(data.size,data.size+data.size*0.1,100)[:,None]
    test_start = svr_rbf.fit(X, data).predict(start)
    test_end = svr_rbf.fit(X, data).predict(end)
    plt.plot(start)
    plt.plot(end)
    plt.show()
    return test_start, test_end

def mirror_extrema_spline(data):
    """
    """

    min_env,max_env,first,last = peak(data)

    if len(min_env) < 2 or len(max_env) < 2:
        finish = True
        top = bot = 0
    else:
        finish = False
        if len(min_env) < 3:
            order_min = 1 # Do linear interpolation if not enough points
        elif len(min_env) < 4:
            order_min = 2 # Do quad interpolation if not enough points
        else:
            order_min = 3

        if len(max_env) < 3:
            order_max = 1  # Do linear interpolation if not enough points
        elif len(max_env) < 4:
            order_max = 2 # Do quad interpolation if not enough points
        else:
            order_max = 3

        start = np.zeros([len(data)/2 - first])
        mid = np.zeros([len(data) - first - (len(data) - last)])
        end = np.zeros([len(data)/2 - (len(data) - last)])

        start = data[first:len(data)/2][::-1]
        mid = data[first:last]
        end = data[len(data)/2:last][::-1]

        mirror = np.append(start,mid)
        mirror = np.append(mirror,end)

        min_env,max_env,first1,last1 = peak(mirror)

        t = interpolate.splrep(min_env, mirror[min_env],
                               k=order_min, per=1)
        top = interpolate.splev(
                            np.arange(len(mirror)), t)

        b = interpolate.splrep(max_env, mirror[max_env],
                               k=order_max, per=1)
        bot = interpolate.splev(
                            np.arange(len(mirror)), b)

        top = top[len(start)-first - len(data) + last:len(start)+len(mid)]
        bot = bot[len(start)-first - len(data) + last:len(start)+len(mid)]

        x = np.arange(len(mirror))
        plt.plot(x,mirror)
        plt.plot(x[len(start)-first - len(data) + last:len(start)+len(mid)],top)
        plt.plot(x[len(start)-first - len(data) + last:len(start)+len(mid)],bot)
        plt.axvline(len(start)-first - len(data)+last)
        plt.axvline(len(start)+len(mid))
        plt.show()

    return top,bot,finish

def emd(data, nimfs=12, shifting_distance=0.2):

        base = len(data)
        signals = np.zeros([base, 2])
        nimfs = range(nimfs) # Max number of IMFs
        IMFs = np.zeros([base, len(nimfs)])
        ncomp = 0
        residual = data
        signals[:, 0] = data

        for j in nimfs:
            # Extract at most nimfs IMFs no more IMFs to be found if Finish is True
            k = 0
            sd = 1.
            finish = False

            while sd > shifting_distance and finish == False:
                print "While j=%i k=%i"%(j,k)
#                import pdb; pdb.set_trace()
                top,bot,finish = mirror_extrema_spline(signals[:, 0])
                print finish
                if finish == False:
                #Calculate the Mean and remove from the data set.
                    mean = (top + bot)/2
                    signals[:,1] = signals[:,0] - mean
                    #Calculate the shifting distance which is a measure of
                    #simulartity to previous IMF
                    if k > 0:
                        sd = (np.sum((np.abs(signals[:,0] - signals[:,1])**2))
                        / (np.sum(signals[:,0]**2)))
                        print sd
                    #Set new iteration as previous and loop
                    signals = signals[:,::-1]
                    k += 1


            if finish:
                #If IMF is a straight line we are done here.
                IMFs[:,j]= residual
                ncomp += 1
                break
            elif not(finish):
                IMFs[:,j] = signals[:,0]
                residual = residual - IMFs[:,j]#For j==0 residual is initially data
                signals[:,0] = residual
                ncomp += 1

        return IMFs[:,0:ncomp]

if __name__=="__main__":
    plt.ion()
    basis = np.linspace(0,5000,1000)
    a = 50. * np.sin(basis/18.)
    b = 100. * np.sin(basis/30.)
    period = (basis[1]-basis[0])/60.
    nose = 100 * np.random.rand(1000)
    data = int_data = area_data = a + b + nose
    data_nose = data + nose
    data_nose -= np.mean(data_nose)
    time = np.arange(0,len(basis)*period,period)
    data -= np.mean(data)

    xpoints, ypoints, cc = peak_para(data)