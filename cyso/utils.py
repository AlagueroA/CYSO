# src/utils.py

import numpy as np
import scipy.constants as sc
from astropy.coordinates import SkyCoord
from cv2 import getRotationMatrix2D,warpAffine,BORDER_CONSTANT,INTER_LANCZOS4

## USEFUL CONSTANTS
sigma_to_FWHM = 2.0 * np.sqrt(2.0 * np.log(2))
FWHM_to_sigma = 1.0 / sigma_to_FWHM
arcsec = np.pi / 648000


## CONVERSIONS
def Wm2_to_Jy(nuFnu, nu):
    '''
    Converts from W.m-2 to Jy .
    
    nuFnu : float
        Flux in [W.m-2/pixel]
        
    nu : float
        Frequency in [Hz]
        
    copied from cpinte/pymcfost/utils.py
    '''
    return 1e26 * nuFnu / nu


def Jy_to_Wm2(nuFnu, nu):
    '''
    Converts from Jy to W.m-2 .

    nuFnu : float
        Flux in [W.m-2/pixel]
        
    nu : float
        Frequency in [Hz]
        
    copied from cpinte/pymcfost/utils.py
    '''
    return 1e-26 * nuFnu * nu



def Wm2_to_ph(nuFnu, nu, surface, exp_time):
    '''
    Converts flux from W.m-2.px-1 to a number of photons .
    
    nuFnu : float
        Flux in [W.m-2/pixel]
        
    nu : float
        Frequency in [Hz]
        
    surface : float
        Telescope surface in [m2]
        
    exp_time : float
        Exposure time of a single frame in [s]
    '''
    E = sc.h*nu
    
    N_ph = nuFnu * exp_time * surface * 1 / E   #multiply by 1 because N_ph in 1px

    return( N_ph )


def ph_to_Wm2(N_ph, nu, surface, exp_time):
    '''
    Converts from a number of photons to W.m-2.px-1 .
    nu : float
        Frequency in [Hz]
        
    surface : float
        Telescope surface in [m2]
        
    exp_time : float
        Exposure time of a single frame in [s]
    '''
    E = sc.h*nu
    
    I = N_ph * E / (exp_time*surface*1)
    return( I )




## OTHERS
def ref_flux(wavelength):
    '''
    Returns the reference flux in W.m-2 to calculate a magnitude in the Johnson system. Filters are taken from ESO Skycalc.

    wavelength : float
        Wavelength in [µm]

    pixel_scale : float
        Pixel scale in [arcsec]

    '''

    #filters taken from ESO Skycalc,  [name, wavelength in nm, flux in erg.s-1.cm-2.A-1]
    filters = [
    ['U',     365, 4.18023e-09],
    ['B',    440, 6.60085e-09],
    ['V',    550, 3.60994e-09],
    ['R',    650, 2.28665e-09],
    ['I',   800, 1.22603e-09],
    ['Z',    900, 7.76068e-10],
    ['Y',    1025, 5.973e-10],
    ['J',    1220, 3.12e-10],
    ['H',     1630, 1.14e-10],
    ['K',     2190,  3.94e-11],
    ['L',     3450, 4.83e-12],
    ['M',     4750, 2.04e-12],
    ['N',     10200, 1.23e-13],
    ['Q',     21000, 6.8e-15]
    ]
    for i in range(len(filters)):                  #flux conversion in W.m-2
        filters[i][2] = 1e-3 * filters[i][1]*10 * filters[i][2]

    #looking for the right filter
    ind_filt = 0
    min_ = filters[0][1]
    for i in range(len(filters)):
        delta = np.abs(wavelength-filters[i][1]/1000)
        if delta < min_ :
            min_ = delta
            ind_filt = i

    #corresponding flux and intensity
    F0 = filters[ind_filt][2]

    print('Reference flux found for band',filters[ind_filt][0])

    return F0


def para_angle(nb_im, exp_time, lat, radec):
    '''
    Returns a list of parallactic angles from the system's coordinates and the exposure time.
    
    nb_im : int
        Number of frames in the sequence
        
    exp_time : float
        Exposure time of a single frame in [s]
        
    lat : float
        Lattitude of the system in [deg]

    radec : SkyCoord object
        Astronomical coordinates of the system

    ##Returns
    pa : list of floats
        List of parallactic angles in [deg]
    '''

    #
    radec_ = SkyCoord(radec, frame='icrs')
    
    #conversion in rad
    lat = lat * np.pi/180
    ra  = radec_.ra.degree * np.pi/180
    dec = radec_.dec.degree * np.pi/180
    
    #time sequence (hour angle), centred on zenith passage at the middle of the sequence
    ha = np.linspace(-exp_time*nb_im//2, exp_time*(nb_im//2+nb_im%2), nb_im)  #in s
    ha = ha / (3600/15) #in h
    ha = ha * np.pi/180 #in rad

    #parallactic angle - find more at https://en.wikipedia.org/wiki/Parallactic_angle
    pa = np.zeros_like(ha)
    for frame in range(nb_im):
        num_tanq = np.sin(ha[frame])
        denum_tanq = np.cos(dec)*np.tan(lat) - np.sin(dec)*np.cos(ha[frame])

        q = np.arctan2(num_tanq, denum_tanq)
        q = q * 180/np.pi

        #q in [-pi; pi]
        if q > 90 :
            q -= 180

        pa[frame] = q

    return(pa)



def prepare_PSF(PSF, FoV_PSF, image, crop=False, npix_crop=None, norm_method='sum', norm_factor=1):
    '''
    Resize, normalize and crop the given PSF sequence to the same pixelscale and size of the given image object
    
    FoV_PSF in ", can be a list-like if rectangular FoV
    
    norm_method : sum or max
    
    crop_factor : length of the original image / length of the cropped image side
    
    '''
    
    try :
        FoV_x = FoV_PSF[0]
        FoV_y = FoV_PSF[1]
    except TypeError :
        FoV_x, FoV_y = FoV_PSF, FoV_PSF
        
    #PSF dimensions
    PSF = np.array(PSF)  #convert PSF into a np array to get the shape attribute
    PSF_list = (len(PSF.shape) == 3)  #checking if a PSF list or a unique PSF is given
    nx_psf, ny_psf = PSF.shape[-2], PSF.shape[-1]
    
    print('Old PSF dimensions : ('+str(nx_psf)+', '+str(ny_psf)+')')
    
    #New PSF dimensions - nb of pixels necessary for the image FOV to be FOV_psf
    new_nx_psf = int(round(FoV_x / image.pixelscale))
    new_ny_psf = int(round(FoV_y / image.pixelscale))
    
    print('New PSF dimensions : ('+str(new_nx_psf)+', '+str(new_ny_psf)+')')
    
    #Resizing
    print('Resizing PSF ...')
    scale = 1e40  #used to avoid approximation errors
    if PSF_list :
        new_PSF = np.zeros((PSF.shape[0], new_nx_psf, new_ny_psf))
        for ind_psf in range(len(PSF)):
            new_PSF[ind_psf] = cv2.resize(PSF[ind_psf]*scale, [new_nx_psf, new_ny_psf], interpolation=cv2.INTER_LINEAR) /scale
        PSF = copy.deepcopy(new_PSF)
    else :
        PSF = cv2.resize(PSF*scale, [new_nx_psf, new_ny_psf], interpolation=cv2.INTER_LINEAR) /scale
    
    
    #Normalizing
    print('Normalizing by '+norm_method)
    if norm_method == 'sum' :
        n_method = np.sum
    elif norm_method == 'max' :
        n_method = np.max
        
    if PSF_list :
        for ind_psf in range(len(PSF)):
            PSF[ind_psf] /= n_method(PSF[ind_psf])
            PSF[ind_psf] *= norm_factor
    else :
        PSF /= n_method(PSF)
        PSF *= norm_factor

            
    #Cropping
    if crop :
        print('Cropping PSF ...')
        
        #get PSF center with max(PSF)
        #new_cx_psf, new_cy_psf = np.unravel_index(PSF.argmax(), PSF.shape)[-2], np.unravel_index(PSF.argmax(), PSF.shape)[-1]
        cx_psf, cy_psf = PSF.shape[-2]//2 + PSF.shape[-2]%2, PSF.shape[-1]//2 + PSF.shape[-1]%2
        print('PSF center is found in coordinate :', cx_psf, cy_psf)
        
        #adjust crop_factor
        offcrop_x, offcrop_y = self.scene_nx//2, self.scene_nx//2  #image is square
        rest_x, rest_y = self.scene_nx%2, self.scene_nx%2
        
        #Centering
        x_min = int(round( cx_psf -  offcrop_x ))
        x_max = int(round( cx_psf +  offcrop_x+rest_x ))
        y_min = int(round( cy_psf -  offcrop_y ))
        y_max = int(round( cy_psf +  offcrop_y+rest_y ))
    
        if x_min < 0 or y_min < 0 :
            print('WARNING : Cropping too harsh, may not work')
        
        #Slicing
        temp = np.zeros((self.image.shape[0], x_max-x_min, y_max-y_min))
        for frame in range(len(PSF)):
            new_PSF[ind_psf] = PSF[ind_psf, x_min:x_max, y_min:y_max]
        else :
            PSF = PSF[x_min:x_max, y_min:y_max]
        PSF = copy.deepcopy(new_PSF)
    
        print('New PSF dimensions after cropping : ('+str(PSF.shape[-1])+', '+str(PSF.shape[-2])+')')
            
    return(PSF)



def find_secondmax(arr):
    '''Return the position of the 2nd maximum in an array'''
    
    flattened = arr.flatten()


    unique_vals = np.unique(flattened)
    sorted_vals = np.sort(unique_vals)

    second_max_val = sorted_vals[-2]

    position = np.argwhere(arr == second_max_val)
    
    return position


def frame_rotate(arr,angle,cx,cy):
    '''rotate a 2D array from an angle using open-cv'''

    M = getRotationMatrix2D((cx, cy), angle, 1)
    array_rot = warpAffine(arr.astype(np.float32), M, (arr.shape[0], arr.shape[1]), flags=INTER_LANCZOS4, borderMode=BORDER_CONSTANT)

    return array_rot
