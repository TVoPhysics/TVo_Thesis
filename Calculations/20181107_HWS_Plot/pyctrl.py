"""
some control/signal processing related functions
"""

import numpy as np
import scipy.signal as sig
import scipy.interpolate as interp
import h5py as h5
from copy import deepcopy

import nds2

### reading data from frame
def fetch_data_nds(gps_start, duration, chan_list, fs=-1, \
                   ifo='H1', data_dir='data/', \
                   server='nds.ligo.caltech.edu', port=31200):
    conn = nds2.connection(server, port)
    chan_list_nds=[]
    for i in range(len(chan_list)):
        chan_list_nds.append('%s:%s'%(ifo, chan_list[i]))
    
    print 'fetching data...'
    data_buffer=conn.fetch(gps_start, gps_start+duration, chan_list_nds)
    print 'fetched!!!'
    
    print 'saving data to: %s%s_%i_%i.h5'%(data_dir, ifo, gps_start, duration)
    f=h5.File('%s%s_%i_%i.h5'%(data_dir, ifo, gps_start, duration), 'a')
    
    for i in range(len(chan_list)):
        
        data=data_buffer[i]
    
        timeSeries=data.data
        fs0=data.sample_rate
        
        print i
        print chan_list[i]
        
        if fs>0:
            if i==0:
                print 'down-sampling to %f Hz'%fs
            fs0 = data.sample_rate
            down_factor = int(fs0/fs)
            fir_aa = sig.firwin(20*down_factor+2, 0.8/down_factor, \
                           window='blackmanharris')
            timeSeries=sig.decimate(timeSeries, down_factor, \
                           ftype=sig.dlti(fir_aa[1:-1], 1), zero_phase=True)
        else:
            if i==0:
                print 'no down-sampling; fs=%f Hz'%fs0
            fs=fs0
        
        if chan_list[i] in f.keys():
            del f[chan_list[i]]     
        grp=f.create_group('%s'%chan_list[i])
        grp.create_dataset('data', timeSeries.shape, dtype='<f8', data=timeSeries)
        grp['t0']=gps_start
        grp['fs']=fs
        
    f.close()
    

### reading time series data from local h5 files
def read_data(gps_start, duration, chan, \
        ifo='H1', data_dir='data/'):
    f=h5.File('%s%s_%i_%i.h5'%(data_dir, ifo, gps_start, duration), 'r')
    grp=f[chan]
    data=grp['data'].value
    fs=grp['fs'].value
    return data, fs

def read_data_multi_chan(gps_start, duration, chan_list, \
        ifo='H1', data_dir='data/'):
    nChan=len(chan_list)
    for i in range(nChan):
        data, fs=read_data(gps_start, duration, chan_list[i], \
                     ifo=ifo, data_dir=data_dir)

        if i==0:
            dataList=np.zeros([nChan, len(data)])

        dataList[i, :]=data
    return dataList, fs

### reading f-domain tf
def read_tf(freq, fName):
    data=np.loadtxt(fName)
    f_in, mag_in, phase_in=data[:, 0], data[:, 1], data[:, 2]
    mag_func=interp.interp1d(np.log(f_in), np.log(mag_in), \
            kind='linear', bounds_error=False, fill_value='extrapolate')
    mag=np.exp(mag_func(np.log(freq)))
    
    phase_func=interp.interp1d(np.log(f_in), phase_in, 
            kind='linear', bounds_error=False, fill_value='extrapolate')
    phase=phase_func(np.log(freq))
    
    tf=mag*np.exp(1j*phase)
    return tf
    
def read_mag(freq, fName):
    data=np.loadtxt(fName)
    f_in, mag_in=data[:, 0], data[:, 1]
    mag_func=interp.interp1d(np.log(f_in), np.log(mag_in), \
            kind='linear', bounds_error=False, fill_value='extrapolate')
            
    mag=np.exp(mag_func(np.log(freq)))
    return mag

## discrete <-> continuous conversion
def d2c(zpk_z, fs, method='tustin', f_match=0.):
    """
    cf. https://www.mathworks.com/help/control/ug/continuous-discrete-conversion-methods.html#bs78nig-12
    """
    zpk=deepcopy(zpk_z)
    zz, pz, kz=zpk
    Ts=1./fs
    
    nzz, npz=len(zz), len(pz)
    zs, ps=np.zeros(nzz, dtype=np.complex), np.zeros(npz, dtype=np.complex)
    ks=kz
    
    if method.lower() == 'tustin':
        zs=(2./Ts)*(zz-1.)/(zz+1.)
        ps=(2./Ts)*(pz-1.)/(pz+1.)
        
        for i in range(nzz):
            ks*=(1.+zz[i])
        for i in range(npz):
            ks/=(1.+pz[i])
            
        if npz>nzz:
            zs_pad=np.ones(npz-nzz, dtype=np.complex)
            zs_pad*=(2./Ts)
            zs=np.hstack([zs, zs_pad])
            ks*=(-1)**(npz-nzz)
        elif nzz>npz:
            ps_pad=np.ones(nzz-npz, dtype=np.complex)
            ps_pad*=(2./Ts)
            ps=np.hstack([ps, ps_pad])
            ks*=(-1)**(nzz-npz)
            
    else: # use direct matching, i.e., the `matched' method in matlab
        zs=np.log(zz)/Ts
        ps=np.log(pz)/Ts
        
        __, k0=sig.freqresp((zs, ps, 1), 2.*np.pi*f_match)
        __, k1=sig.freqz_zpk(zz, pz, kz, np.pi*f_match/(fs/2.))
        ks=k1/k0
        
    return (zs, ps, ks)

def c2d(zpk_s, fs, method='tustin', f_match=0):
    zpk=deepcopy(zpk_s)
    zs, ps, ks=zpk
    Ts=1./fs
    
    nzs, nps=len(zs), len(ps)
    zz, pz=np.zeros(nzs, dtype=np.complex), np.zeros(nps, dtype=np.complex)
    
    kz=ks
    
    if method.lower() == 'tustin':
        zz=(1.+0.5*Ts*zs)/(1.-0.5*Ts*zs)
        pz=(1.+0.5*Ts*ps)/(1.-0.5*Ts*ps)
        
        for i in range(nzs):
            kz*=(1.-0.5*Ts*zs[i])
        for j in range(nps):
            kz/=(1.-0.5*Ts*ps[j])
            
        kz*=((2./Ts)**float(nzs-nps))
        
        if nps>nzs:
            zz_pad=-np.ones(nps-nzs, dtype=np.complex)
            zz=np.hstack([zz, zz_pad])
        elif nzs>nps:
            pz_pad=-np.ones(nzs-nps, dtype=np.complex)
            pz=np.hstack([pz, pz_pad])
            
    else:
        zz=np.exp(zs*Ts)
        pz=np.exp(ps*Ts)
        
        __, k0=sig.freqz_zpk(zz, pz, 1, 2.*np.pi*f_match)
        __, k1=sig.freqresp(zpk, np.pi*f_match/(fs/2.))
        kz=k1/k0
        
    return (zz, pz, kz)
    
    
# discard similar zpk pairs
def discard_similar_pz(zpk_s, cut=0.05, f_match=0.):
    zz, pp, kk=zpk_s
    
    z1=np.where(np.imag(zz)>0.)
    z1=zz[z1]
    
    z2=np.where(np.imag(zz) == 0.)
    z2=zz[z2]
    
    p1=np.where(np.imag(pp)>0.)
    p1=pp[p1]
    
    p2=np.where(np.imag(pp) == 0.)
    p2=pp[p2]
    
    np1=len(p1)
    nz1=len(z1)
    
    z_del=[]
    p_del=[]
    for i in range(nz1):
        j=np.argmin(np.abs(np.abs(p1)-np.abs(z1[i])))
        if np.abs(np.abs(p1[j])-np.abs(z1[i]))<cut and not j in p_del:
            z_del.append(i)
            p_del.append(j)
            
            kk*=(f_match-z1[i])*(f_match-np.conj(z1[i]))\
                /(f_match-p1[j])/(f_match-np.conj(p1[j]))
            
    z1=np.delete(z1, z_del)
    p1=np.delete(p1, p_del)
    
    idx_z=np.argsort(np.abs(z1))
    z1 = z1[idx_z]
    idx_p=np.argsort(np.abs(p1))
    p1 = p1[idx_p]
    
    zz = np.zeros(2 * len(z1) + len(z2), dtype=np.complex)
    pp = np.zeros(2 * len(p1) + len(z2), dtype=np.complex)
    
    for i in range (len(z1)):
        zz[2*i] = z1[i]
        zz[2*i+1] = np.conj(z1[i])
        
    for i in range(len(p1)):
        pp[2*i] = p1[i]
        pp[2*i+1] = np.conj(p1[i])
        
    zz[2*len(z1):]=z2
    pp[2*len(p1):]=p2

    return (zz, pp, kk)
    
# res_g
def get_res_g_pole_pair(f0, Q):
    w0=2.*np.pi*f0
    a=-w0/2./Q
    b=np.sqrt(w0**2.-a**2.)
    
    a=np.round(a, 6)
    b=np.round(b, 6)
    
    return a, b

def get_res_f0_Q(a, b):
    w0=np.sqrt(a**2.+b**2.)
    f0=w0/2./np.pi
    Q=-w0/2./a
    return f0, Q
    

# foton
def get_foton(zpk_s):
    """
    print out a zpk_s=(zz, pp, kk) (in s-domain) in the foton format
    """
    zz, pp, kk=zpk_s
    print('zpk(')
    
    print('[')
    for z in zz:
        if np.imag(z)>=0:
            print('%e+i*%e; '%(np.real(z),np.imag(z)))
        else:
            print('%e-i*%e; '%(np.real(z),-np.imag(z)))
    print('],')
    
    print('[')
    for p in pp:
        if np.imag(p)>=0:
            print('%e+i*%e; '%(np.real(p),np.imag(p)))
        else:
            print('%e-i*%e; '%(np.real(p),-np.imag(p)))
    print('],')
    
    print('%e'%np.real(kk))
    print(')')
        

# MISO coherence and TF
def mcoherence(tar, wit, \
               fs=1., Tperseg=1.):
    """
    input: 
        tar = [nTime]; an 1d vector for the target timeseries
        wit = [nChan, nTime]; first index for channel, second for time stamp; timeseries for the witness channels
        fs = scalar; sampling frequency in [Hz]
        Tperseg = scalar; fft length in [sec]
    output:
        freq = [nFreq]; 1d frequency vector
        coh = [nFreq]; 1d vector; multi-witness coherence evaluated at each freq
        tf = [nChan, nFreq]; transfer functions. 
    """
    nperseg = int(fs*Tperseg)
    if len(wit.shape)==1:
        nChan = 1; nTime=len(wit)
        wit=np.reshape(wit, (nChan, nTime))
    elif len(wit.shape)==2:
        nChan, nTime=wit.shape
    else:
        raise ValueError('Input wit got wrong shape! It should have [nChan, nTime]!')
    
    freq, Pbb = sig.csd(tar, tar, fs=fs, nperseg=nperseg)
    nFreq=len(freq)
    
    tf = np.zeros([nChan, nFreq], dtype=np.complex)
    Paa = np.zeros([nChan, nChan, nFreq], dtype=np.complex)
    Pab = np.zeros([nChan, nFreq], dtype=np.complex)
    coh =np.zeros(nFreq)
            
    Paa_inv = np.zeros([nChan, nChan, nFreq], dtype=np.complex)
    
    for i in range(nChan):
        freq, Pab_single = sig.csd(tar, wit[i, :], fs=fs, nperseg=nperseg)
        
        if i==0:
            nFreq=len(freq)
            
        Pab[i, :]=Pab_single
        
        for j in range(i, nChan, 1):
            freq, Paa_single = sig.csd(wit[i, :], wit[j, :], fs=fs, nperseg=nperseg)
            Paa[i, j, :] = Paa_single
            Paa[j, i, :] = np.conj(Paa_single)
            
    for k in range(nFreq):
        Paa_inv[:, :, k] = np.linalg.inv(Paa[:, :, k])
        tf[:, k] = np.dot(Pab[:, k], Paa_inv[:, :, k])
        #coh[k] = np.abs(np.dot(np.conj(tf[:, k]), Pab[:, k]))/Pbb[k]
        coh[k] = np.real(np.dot(tf[:, k], np.conj(Pab[:, k]))/Pbb[k])
        
    return freq, coh, tf

def get_spec(data, \
        fs=64., Tperseg=64, Ttrunc=0):
    nperseg=int(fs*Tperseg)
    ntrunc=int(fs*Ttrunc)
    freq, P=sig.welch(data[ntrunc:], fs=fs, nperseg=nperseg, axis=-1)
    rms=rms_vs_freq(P, freq)
    return freq, P, rms

def rms_vs_freq(P, freq):
    df=freq[1]-freq[0]
    rms=np.zeros(P.shape)
    try:
        nChan, nFreq=P.shape
        for i in range(nFreq):
            rms[:, i]=np.sqrt(np.sum(P[:, i:], axis=-1)*df)
    except ValueError:
        nFreq=len(freq)
        for i in range(nFreq):
            rms[i]=np.sqrt(np.sum(P[i:])*df)
    return rms

            
