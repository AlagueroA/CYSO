# src/CYSO_input.py

## IMAGE
path_image = ''
mcfost = False
radmc3d = False

wavelength = None   #wavelength in [µm], put None if you use mcfost or radmc3d
pixelscale = None   #pixel scale in [arcsec/pixel], put None if you use mcfost or radmc3d
unit = None         #brightness unit, put None if you use mcfost or radmc3d


## SYSTEM
distance = 140   #pc

path_pa = ''     #path to parallactic angles fits file, in [xxx], let empty to create one
radec = '14h08m10.155s -41d23m52.57s'  #ra dec of the object in hhmmss and ddmmss in the icrs frame, used if you let path_pa empty
lat = -24.6280555          #lattitude of the observer in [deg], used if you let path_pa empty



## INSTRUMENT
coro = 20             #coronagraph size in [mas], put to 0 if none
tel_diam = 39         #telescope diameter in [m]
tel_surface = 1199    #telescope surface in [m2]

phot_noise = True
readout_noise = 20    #in [photons], put to 0 if none
atm_bg_noise = 15     #in [mag/arcsec2], put to 0 if none



## PSF SEQUENCE
path_PSF_onaxis = ''    #let empty if None
path_PSF_offaxis = ''   #let empty if None
nb_frames = 1   #number of frames to use
exp_time = 10   #exposure time of each individual frames in [s]


combine_method = 'median'    #recombination method

pixelscale_PSF = 6/512  #pixel scale in [arcsec/pixel]
FoV = 6  #field of view in [arcsec]
norm_cst = 1   #normalisation constant for the PSF

## POST-PROCESSING
RDI = False
path_image_RDI = ''
path_PSF_onaxis_RDI = ''   #let empty to use same as path_PSF_onaxis
path_PSF_offaxis_RDI = ''   #let empty to use same as path_PSF_offaxis

ADI = False


## EXPORT
export_fits = True
