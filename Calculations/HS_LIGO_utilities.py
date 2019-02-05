"""
HWS Utility Tools

This is a collection of tools for interacting with the Hartmann wavefront sensors.
This mostly functions for getting and plotting data.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import h5py
import glob
import shutil
import sys
import os
import os.path
import subprocess
import pickle
import ezca

from HWS.HS_Image import HS_Image
from HWS.HS_Centroids import HS_Centroids
from HWS.HS_Gradients import HS_Gradients
from HWS.HS_WFP import HS_WFP
from HWS.HS_Export import HS_Export

import HWS.HSM_WFN

from collections import OrderedDict

ez = ezca.Ezca()

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def caget(channel):
    rtn = subprocess.check_output(['caget', '-g', '15', channel])
    return rtn.split()[-1]

def getFileGPSTimes(fList):
    """
    For a list of HWS file names this extracts the GPS times
    from the file names and returns them as a list
    
    Parameters
    ==========
    flist : list
        List of paths to HWS files
    
    Returns
    =======
    GPS times, np.array[int]
    """
    gpstimeList = []
    for f in fList:
        gpstimeList.append(int(os.path.basename(f).split('_')[0]))
    return np.array(gpstimeList, dtype=np.int)


def date_2_gps(date):
    """
    Calls the gpstime command line tool to convert a date string
    into GPS time

    Returns
    =======
    GPS time, float

    Example
    =======
    date_2_gps('2018-11-01 17:36:55.549957 UTC')
    # Returns 1225129033.55
    """
    rtn = subprocess.check_output(['gpstime', date])
    return float(rtn.split("GPS:")[-1].strip())


def gps_2_date(date):
    """
    Calls the gpstime command line tool to convert a gps time
    into UTC date

    Returns
    =======
    Date, datestring
    """
    rtn = subprocess.check_output(['gpstime', date])
    return float(rtn.split("GPS:")[-1].strip())

def any_2_gps(date_or_gps):
    """
    Always returns a GPS time. If a date is entered then
    date_2_gps is called to convert it.
    
    Returns
    =======
    GPS time, float
    """
    try:
        gps = float(date_or_gps)
    except ValueError:
        gps = date_2_gps(date_or_gps)
    
    return gps


def _look_for_data_dir():
    """
    Looks for some default directories that are L1 and H1 specific
    where HWS data is stored.

    Returns
    =======
    Path to HWS data : str
    """

    ifo = os.environ['IFO']

    if ifo == 'H1':
        if os.path.exists('/h1hwsmsr/data'):
            data_dir = '/h1hwsmsr/data'
        elif os.path.exists('/data'):
            data_dir = '/data'
        else:
            raise FileNotFoundError("Could not find default data directory")
    elif ifo == 'L1':
        raise Exception("TODO, check LLO storage locations")
    elif ifo == 'TEST':
        data_dir = '/home/ddb/hws/data'
    else:
        raise Exception("Unexpected IFO input `{}`, should be H1 or L1".format(ifo))

    return data_dir


def get_servers(ifo=None):
    """
    This function returns the HWS servers that are running at
    the site. By default it will return the servers for the
    current site you run it at.

    Parameters
    ----------
    ifo : str, optional
        Interferometer name, either H1 or L1. When set to None it
        tried to read the $IFO environment variable
   
    Returns
    -------
    Dict of servers: dict(optic:(server, port))
    """
    if ifo is None:
        if 'IFO' not in os.environ:
            raise Exception("IFO environment variable is not set")
        else:
            ifo = os.getenv('IFO')


    if ifo == 'H1':
        return    {
                    'ITMX': ('h1hwsmsr', 9991),
                    'ITMY': ('h1hwsmsr1', 9991),
                    'ETMX': ('h1hwsex', 9991),
                    'ETMY': ('h1hwsey', 9991),
                  }
    elif ifo == 'L1':
        # Are these right for L1????
        return    {
                    'ITMX': ('l1hwsmsr', 9991),
                    'ITMY': ('l1hwsmsr1', 9991),
                    'ETMX': ('l1hwsex', 9991),
                    'ETMY': ('l1hwsey', 9991),
                  }
    else:
        raise Exception("Argument {} is not L1 or H1".format(ifo))


def get_file(time, dtime=60, optic='ITMX', data_dir=None, ifo=None,  extension='.hdf5'):
    """
    Searchs for a HWS file at a certain GPS time or date.

    Parameters
    ==========
    time : int or str
        GPS time or date to search for
    dtime : int, optional
        Bound to search within in seconds from time
    data_dir : str, optional
        Top level data storage directory
    ifo : str, optional
        Interferometer, H1 or L1
    optic : str, optional
        The choice of optic to look at, ITMX, ETMX, ITMY, ETMY
    extension : str, optional
        Type of file to look for
        
    Returns
    =======
    GPS time, file names, tuple
    """
    if ifo is None:
        ifo = os.environ['IFO']
    
    assert(ifo in ['H1','L1','TEST'])
    
    if data_dir is None:
        data_dir = _look_for_data_dir()
    data_dir = os.path.join(data_dir, ifo, optic+"_HWS/")
    t0       = any_2_gps(time)
    fileList = []
    
    print(data_dir)
    
    # Search in nearest time folders just in case 
    for i in [-1,0,1]:
        tSTEM = str(int(np.floor(t0/1E5)) + i) + '00000'
        fileTemplate = os.path.join(data_dir, tSTEM, '*' + extension)
        fileList.extend( glob.glob(fileTemplate) )

    if len(fileList) == 0: raise FileNotFoundError("No files for time {} +- {} seconds".format(time, dtime))
    fileList = np.array(fileList)
    
    # find the closest file to input file time
    fileTimes   = np.array( getFileGPSTimes(fileList) )
    idx = abs(fileTimes-time).argmin()
    
    return fileTimes[idx], fileList[idx]


def get_files(start, end, optic='ITMX', data_dir=None, ifo=None, extension='.hdf5', dtime=60):
    """
    Gets all the HWS files between two GPS times. Can be customised
    to run at different sites and storage locations

    Parameters
    ==========
    start : int or str
        Start GPS time or date
    end : int or str
        End GPS time or date
    optic : str, optional
        The choice of optic to look at, ITMX, ETMX, ITMY, ETMY
    ifo : str, optional
        Interferometer, H1 or L1. If none provided it will use the $IFO environment variable
    data_dir : str, optional
        Top level data storage directory. If None then it will guess some likely options:
            H1, 1) /h1hwsmsr/data
            H1, 2) /data
            L1, 1) ???? 
    dtime : int, optional
        Search within this time bound of start and end for a file
        
    Returns
    =======
    Dictionary of GPS times and file names, key: GPS times, value: filename
    """
    if ifo is None:
        ifo = os.environ['IFO']
    
    assert(ifo in ['H1','L1'])
    
    if data_dir is None:
        data_dir = _look_for_data_dir()
        
    HWStimes = OrderedDict()
    data_dir = os.path.join(data_dir, ifo, optic+"_HWS/")
    start_gps = any_2_gps(start)
    end_gps   = any_2_gps(end)
    
    # We need to get all files between two times, so we may have to check multiple folders
    start_folder = int(np.round(start_gps/1e5))
    end_folder   = int(np.round(end_gps/1e5)  )
    dfold = end_folder - start_folder
    
    fileList = []
    
    for i in list(range(int(start_folder)-1, int(end_folder)+2)):
        tSTEM = str(i*100000)
        fileTemplate = os.path.join(data_dir, tSTEM, '*' + extension)
        fileList.extend( glob.glob(fileTemplate) )

    fileList = np.array(fileList)
    if len(fileList) == 0:
        raise FileNotFoundError("No HWS files found between {} and {} (ext={})".format(int(start), int(end), extension))

    # find the closest file to input file time
    fileTimes = np.array( getFileGPSTimes(fileList) )
    idx       = np.array(np.logical_and(fileTimes>=(start_gps-dtime), fileTimes<=(end_gps+dtime)))
    
    if len(fileTimes[idx]) == 0:
        raise FileNotFoundError("No HWS files found between {} and {} (ext={})".format(int(start), int(end), extension))
        
    return OrderedDict(sorted(zip(fileTimes[idx], fileList[idx])))


def copy_files(destination, start, end, *args, **kwargs):
    """
    Copies HWS files from a data directory to some storage location
    
    Parameters
    ==========
    destination : str
        Path to copy files to
    start : int or str
        Start GPS time or date
    end : int or str
        End GPS time or date
    data_dir : str
        Top level data storage directory
    ifo : str
        Interferometer, H1 or L1
    optic : str
        The choice of optic to look at, ITMX, ETMX, ITMY, ETMY
    dtime : int, optional
        Search within this time bound of start and end times for files
        
    Returns
    =======
    Dict(GPS times: copied files)
    """
    if not os.path.exists(destination):
        os.mkdir(destination)
        
    files = get_hws_files(start, end, *args, **kwargs)
    
    new_files = []
    
    for f in files.values():
        shutil.copy(f, destination)
        new_files.append( os.path.join(destination, os.path.basename(f)) )
    
    return OrderedDict(sorted(zip(files.keys(), new_files)))


def open_file(filename):
    """
    Opens a HWS file. Handles different extensions depending on storage format.
    This is either a pickle file for old storage or HDF5 for the new storage
    format. Either way this returns a `HS_Export` object to be used.

    Parameters
    ==========
    filename : str
        Path and filename to a HWS file

    Returns
    =======
    Loaded data : HS_Export
    """
    ext = os.path.splitext(filename)[1]
    
    if ext == '.p':
        # Stupid hack to get around the fact HWS is now packaged
        # instead of separate scripts
        remove = False
        if 'HS_Export' not in sys.modules:
            remove = True
            sys.modules['HS_Export'] = HWS.HS_Export
        data = pickle.load(open(filename, 'rb'))
        if remove: del sys.modules['HS_Export']
        return data 
    elif ext == '.hdf5':
        
        data = h5py.File(filename, 'r')
        
        exp = HS_Export()
        exp.avGradients = data['avGradients'].value
        exp.avImage = data['avImage'].value

        for attr in data['avGradients'].attrs.keys():
            setattr(exp, attr, data['avGradients'].attrs[attr])

        return exp
    else:
        raise Exception("Unhandled extension")


def get_bad_pixel_file(optic):
    data_dir = _look_for_data_dir()
    ifo = os.environ['IFO']
    fname = os.path.join(data_dir, ifo, "BAD_PIXELS", optic + ".npy")
    if os.path.exists(fname):
        return fname
    else:
        return None


def get_centroids(image, ref_centroids=None, bad_pixels=None, background_level=50,
                    mask_center=None, mask_radius=None):
    """
    Takes an image and processes it to extract the centroids in it.
    This is the first step in computing any Hartman analysis of an image.
    If the image being referenced to another then the `ref_centroids` object
    should be provided so that the same number of centroids are computed.
    
    This code is based on the processing that is done in the `State_3.py` 
    of the LIGO HWS analysis.

    This code does all the image processing such as removing bad pixels
    and setting circular masks if specified.
    
    Parameters
    ==========
    image : np.array
        Camera image to process
    ref_centroids : HS_Centroids
        Reference centroids to use for extracting centroids from `image`
    bad_pixels : np.array[array]
        This is a special bad pixel array object that contains an element]
        called `coords` which lists each bad pixel location.
    background_level : int
        Sets the background level of the image for processing centroids.
        50 is the default level set in `State_3.py` in the HWS code.
    mask_center : np.array[array, ndim=1, size=2]
        Set center of the mask
    mask_radius : int
        Mask radius in pixels

    Returns
    =======
    Centroids : HS_Centroids
    """
    hsi = HS_Image()
    hsi.original_image = image.astype(np.int)
    hsi.background = background_level

    if bad_pixels is not None:
        bad_pix_array = np.load(bad_pixels)
        hsi.bad_pixels = bad_pix_array.item(0)['coords']
        hsi.to_fix_bad_pixels=True

    if mask_radius is not None and mask_center is not None:
        hsi.to_mask_image = True
        hsi.mask_center = mask_center
        hsi.mask_radius = mask_radius
    
    hsi.process_image()

    hsc = HS_Centroids()
    hsc.hsimage = hsi
    
    if ref_centroids is None:
        hsc.find_centroids_from_image()
    else:
        hsc.find_centroids_using_template(ref_centroids.centroids)
    
    return hsc


def get_gradients(ref_centroids, centroids, origin=[511.5, 511.5], magnification=17.5, lever_arm=0.01):
    """
    Compares two sets of centroids to compute the gradients
    between them.

    This

    Parameters
    ==========
    ref_centroids : HS_Centroids
        Reference centroids
    centroids : HS_Centroids
        Centroids being analysed
    origin : np.array[float, shape=(2,)]
        The origin of the system relative to the (0,0) pixel of
        the camera sensor.
    magnification : float
        Magnification of the HWS optical system, e.g. `caget H1:TCS-ETMX_HWS_MAGNIFICATION`
    lever_arm : float
        Lever arm of the HWS system, e.g. `caget H1:TCS-ETMX_HWS_LEVER_ARM`
    
    Returns
    =======
    Gradients : HS_Gradients
    """
    hsg = HS_Gradients(ref_centroids, centroids)
    hsg.lever_arm = lever_arm
    hsg.magnification = magnification
    hsg.origin = np.array(origin, dtype=np.float64)
    hsg.construct_gradients()
    return hsg


def get_seidel(gradients):
    """
    Given a set of gradients this will compute the
    Seidel coefficients.

    Parameters
    ==========
    gradients : HS_Gradients

    Returns
    =======
    Seidel Coefficients : dict
    """
    hswfp = HS_WFP(gradients, order=6)
    hswfp.compute_poly_coeffs()
    hswfp.compute_seidel_coeffs()    
    return hswfp.seidel_coeffs
    
    
def get_optical_depth(start_gradients, final_gradients, magnification=17.5, N=100):
    """
    Computes the wavefront aberration that ocurred between two sets
    of Hartman gradients. It returns the wavefront aberration in
    terms of optical depth in nanometers.
    
    Parameters
    ==========
    start_gradients : HS_Gradients
        Gradients at some initial reference time
    final_gradients : HS_Gradients
        Gradients at a final time
    magnification : float
        Magnification of the HWS optical system, e.g. `caget H1:TCS-ETMX_HWS_MAGNIFICATION`
    N : int
        Generate an NxN wavefront aberration array
        
    Returns
    =======
    Optical depth in nanometers : np.array[float64]
    x, actual units on mirror : np.array[float]
    y, actual units on mirror : np.array[float]
    """
    hrg = start_gradients.gradients
    hsg = final_gradients.gradients
    
    pixel_size = 12e-6
    pixel_number = 1024

    rng = pixel_size * pixel_number * magnification/2
    
    limits = {'x': [-rng, rng],
              'y': [-rng, rng]}
    spc = 2*rng/N
    x_g, y_g = HWS.HSM_WFN.setup_grid(limits, spc)

    U = hsg[:,0] - hrg[:,0]
    V = hsg[:,1] - hrg[:,1]
    
    U = U - np.mean(U)   # remove tilt
    V = V - np.mean(V)
    
    cents = hsg[:,2:4]
    grads = np.transpose( np.array([U,V]) )

    # calculate the wavefront
    wf = -1E9 * HWS.HSM_WFN.calculate_wf(cents, grads, limits, spc)
    return wf, x_g, y_g


def get_origin(optic):
    return (ez['TCS-{}_HWS_ORIGIN_X'.format(optic)], ez['TCS-{}_HWS_ORIGIN_Y'.format(optic)])

def get_mask(optic):
    return {
            'center_x': ez['TCS-{}_HWS_MASK_CENTER_X'.format(optic)],
            'center_y': ez['TCS-{}_HWS_MASK_CENTER_Y'.format(optic)],
            'radius':   ez['TCS-{}_HWS_MASK_RADIUS'.format(optic)]
        }





