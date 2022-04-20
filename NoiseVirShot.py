import numpy as np
import tqdm
import matplotlib.pyplot as plt

def correlate_chunk(data, virTr=1, dt=None):

    ns,nx = data.shape 
    
    if dt==None:
        print('dt is not set. dt assumed to be 0.008s')
        dt=0.008
    
    
    cc_data=np.zeros((2*ns-1,nx))
    
    for ix in range(nx):
        cc_data[:,ix]=np.correlate(data[:,ix], data[:,virTr], mode='full')
        
    
    return cc_data


def correlate_allchunks(data, virTr=1, dt=None, chunk_size=3):

    ns,nx = data.shape
    
    if dt==None:
        print('dt is not set. dt assumed to be 0.008s')
        dt=0.008
        
    t=np.arange(0,ns*dt,dt)
    t_max=t[-1]
    
    print('total record length [s]:',t_max)
    print('chunck size [s]: ', chunk_size)
    
    nchunks=int(t_max/chunk_size)
    ch_ns=int(chunk_size//dt)
    
    print('number of chuncks:', nchunks)
    print('samples per chunk:', ch_ns)
    
    data_chunks=np.zeros((nchunks, ch_ns,nx))
    
    # counters for creating chunks
    it1=0
    it2=ch_ns
    for ich in range(nchunks):
        data_chunks[ich,:,:] = data[it1:it2,:]
        it1=it2+1
        it2=it1+ch_ns
        
    corr_chunks=np.zeros((nchunks, 2*ch_ns-1,nx))
    
    print('Correlating chunks..')
    for ich in tqdm.tqdm(range(nchunks)):
        corr_chunks[ich,:,:] = correlate_chunk(data_chunks[ich,:,:], virTr=virTr, dt=dt)        
    return corr_chunks
        
    
        

def stack_chunks(corr_chunks, dt=0.008, method='Linear'):

    if method=='Linear':
        stack=np.sum(corr_chunks, axis=0)

    # TO DO: add other stacking methods..

    return stack


def plot_stack(data, dt=None, clim=1e6, **kwargs):

    
    nt,nx = data.shape
    
    if dt == None:
        print('dt is not set. dt assumed to be 0.008s')
        dt = 0.008
    if dt != None:
        print('dt =', dt, 's')
    
    
    tmax = (nt*dt)/2
    tmin = -tmax
            
    if kwargs:
        dx = kwargs['dx']
        xmin = 0
        xmax = dx*nx
        xlabel='Distance [m]'
        print('dx=', dx, 'm')

    else:
        print('dx is not given!')
        xmin = 0
        xmax = nx
        xlabel='Trace number'
    
    
    
    plt.figure()
    plt.imshow(data, extent=[xmin, xmax,tmax, tmin], cmap='gray', aspect='auto')
    plt.clim(-clim, clim)
    plt.ylabel('Lag time [s]')
    plt.xlabel(xlabel)

def NoiseVirShot(data, virTr=1, dt=None, chunk_size=3, stacking_method='Linear'):
    '''
    correlated_stack = NoiseVirShot(data, virTr=1, dt=None, chunk_size=3, stacking_method='Linear')

    Input:
        data: input record of seismic noise

        VirTr: trace to be used as a source

        dt: sampling interval in s. If not given, it is assumed that dt=0.008

        chunk_size: window length to be used for the interferometry in seconds, default value is 3.
        
        stacking_method: method to stack chunks, method is Linear as defualt. Other methods to be added later.

    Output:

        Correlated_stack: result of x-correlation
    '''


    
    chunks = correlate_allchunks(data,virTr, dt=dt, chunk_size=chunk_size)
    stack  = stack_chunks(chunks, method=stacking_method)
    
    return stack