# src/synth_obs.py
from astropy.io import fits
from astropy.convolution import convolve_fft
import os
import copy
import random
from vip_hci.preproc import frame_rotate

from .utils import *



class SynthObs:

    def __init__(self,Image, PSF_onaxis=None, PSF_offaxis=None,
                 coro=0,
                 tel_diam=8,
                 tel_surface=163,
                 phot_noise=True,
                 exp_time=10,
                 readout_noise=0,
                 atm_bg_noise=0,
                 combine_method='median',
                 RDI=False,
                 Image_RDI=None,
                 PSF_onaxisRDI=None,
                 PSF_offaxisRDI=None,
                 ADI=False,
                 path_pa=None,
                 radec = SkyCoord('14h08m10.155s -41d23m52.57s', frame='icrs'),
                 lat=-24.6280555,
                 export_fits=False,
    ):
        '''Ah bah là il va falloir décrire toutes ces variables...'''
    
        #Image
        self.Image = Image
        self.nx = Image.nx
        self.ny = Image.ny
        self.wl = Image.wl
        self.freq = Image.freq
        self.pixelscale = Image.pixelscale
        
        #PSFs
        self.PSF_onaxis = PSF_onaxis
        self.PSF_offaxis = PSF_offaxis
        self.nb_frames = PSF_onaxis.nb_frames  #should be the same as off-axis

        self.PSF_onaxisRDI = PSF_onaxisRDI
        self.PSF_offaxisRDI = PSF_offaxisRDI
        
        #Variables
        self.coro = coro
        self.tel_diam = tel_diam
        self.tel_surface = tel_surface
        self.exp_time = exp_time
        
        self.phot_noise = phot_noise
        self.readout_noise = readout_noise
        self.atm_bg_noise = atm_bg_noise
        
        self.combine_method = combine_method
        
        self.RDI = RDI
        self.Image_RDI = Image_RDI
        
        self.ADI = ADI
        self.path_pa = path_pa
        self.lat = lat
        self.radec = radec
        
        #Products
        self.sequence = np.zeros((self.nb_frames, self.nx, self.ny))
        self.synthobs = np.zeros((self.nx, self.ny))    #recombined sequence
        if self.RDI:
            self.sequence_RDI = np.zeros((self.nb_frames, self.nx, self.ny))
        self.export_fits = export_fits

        #Do the things
        self.get_pa()
        self.create_synth_obs()
        if self.export_fits:
            self.export()

    #
    def get_pa(self):
        '''Get the list of parallactic angles. Either from an external file of calculating it.'''
        if self.ADI:
            if self.path_pa == '':
                print('Calculating the list of parallactic angles for an object with (RA,dec):', self.radec, 'observed at a latitude of', self.lat, 'deg')
                self.list_pa = para_angle(nb_im=self.nb_frames, exp_time=self.exp_time, lat=self.lat, radec=self.radec)
            else:
                print('Opening the list of parallactic angles at:', path_pa)
                try:
                    with fits.open(self.path_pa) as hdul:
                        self.list_pa = hdul[0].data
                except OSError:
                    print(f"Cannot open FITS file: {self.file}")
                except Exception as e:
                    print(f"Error reading {self.file}: {e}")
        else:
            self.list_pa = None

    #
    def create_synth_obs(self):
        '''Wrapper for the whole process'''
        
        #sky background
        if self.atm_bg_noise > 0:
            self.Image.image = apply_atm_bg(I=self.Image.image,
                         wavelength=self.wl,
                         pixelscale=self.pixelscale,
                         nx=self.nx,
                         ny=self.ny,
                         sky_mag=self.atm_bg_noise)
            if self.RDI:
                self.Image_RDI.image = apply_atm_bg(I=self.Image_RDI.image,
                             wavelength=self.wl,
                             pixelscale=self.pixelscale,
                             nx=self.nx,
                             ny=self.ny,
                             sky_mag=self.atm_bg_noise)
        
        

        #convolution
        self.sequence = convolution(I=self.Image.image,
                    PSF_on=self.PSF_onaxis.image,
                    PSF_off=self.PSF_offaxis.image,
                    freq=self.freq,
                    telescope_surface=self.tel_surface,
                    exp_time=self.exp_time,
                    RON=self.readout_noise,
                    nb_frames = self.nb_frames,
                    nx = self.nx,
                    ny = self.ny,
                    coro = self.coro/1000, #in arcsec
                    pixelscale = self.pixelscale,
                    ADI = self.ADI,
                    angle_list = self.list_pa)
        if self.RDI:
            self.sequence_RDI = convolution(I=self.Image_RDI.image,
                        PSF_on=self.PSF_onaxis.image,
                        PSF_off=self.PSF_offaxis.image,
                        freq=self.freq,
                        telescope_surface=self.tel_surface,
                        exp_time=self.exp_time,
                        RON=self.readout_noise,
                        nb_frames = self.nb_frames,
                        nx = self.nx,
                        ny = self.ny,
                        coro = self.coro/1000, #in arcsec
                        pixelscale = self.pixelscale,
                        ADI = self.ADI,
                        angle_list = self.list_pa)

        #photon noise
        if self.phot_noise:
            self.sequence = apply_photon_noise(I_list=self.sequence,
                               freq=self.freq,
                               telescope_surface=self.tel_surface,
                               exp_time=self.exp_time,
                               nb_frames=self.nb_frames,
                               nx=self.nx,
                               ny=self.ny)
            if self.RDI:
                self.sequence_RDI = apply_photon_noise(I_list=self.sequence_RDI,
                                   freq=self.freq,
                                   telescope_surface=self.tel_surface,
                                   exp_time=self.exp_time,
                                   nb_frames=self.nb_frames,
                                   nx=self.nx,
                                   ny=self.ny)

        #read-out noise
        if self.readout_noise > 0:
            self.sequence = apply_readout_noise(I_list=self.sequence,
                                freq=self.freq,
                                telescope_surface=self.tel_surface,
                                exp_time=self.exp_time,
                                RON=self.readout_noise,
                                nb_frames=self.nb_frames,
                                nx=self.nx,
                                ny=self.ny)
            if self.RDI:
                self.sequence_RDI = apply_readout_noise(I_list=self.sequence_RDI,
                                    freq=self.freq,
                                    telescope_surface=self.tel_surface,
                                    exp_time=self.exp_time,
                                    RON=self.readout_noise,
                                    nb_frames=self.nb_frames,
                                    nx=self.nx,
                                    ny=self.ny)

        #RDI
        if self.RDI:
            self.sequence = RDI(I_list=self.sequence,
                I_list_RDI=self.sequence_RDI,
                nb_frames=self.nb_frames,
                nx=self.nx,
                ny=self.ny)

        #ADI
        if self.ADI:
            self.sequence = ADI(I_list=self.sequence,
                angle_list=self.list_pa,
                nb_frames=self.nb_frames,
                nx=self.nx,
                ny=self.ny)

        #collapse the cube
        if self.combine_method is not None:
            self.synthobs = recombine(I_list=self.sequence,
                      mode=self.combine_method,
                      nb_frames=self.nb_frames)


    def export(self):
        '''Export the synthetic observation in a FITS file'''

        #create fits
        hdr = fits.Header()
        hdr['EXT0'] = 'Intensity'
        hdr['IUNIT'] = 'W.m-2'

        hdr['NX'] = self.nx
        hdr['NY'] = self.ny
        hdr['WL'] = self.wl
        hdr['WLUNIT'] = 'mu m'
        hdr['FREQ'] = self.freq
        hdr['FREQUNIT'] = 'Hz'

        hdr['NB_FRAMES'] = self.nb_frames
        hdr['IWA'] = self.coro
        hdr['IWAUNIT'] = 'mas'
        hdr['TEL_DIAM'] = self.tel_diam
        hdr['TEL_DIAMUNIT'] = 'm'
        hdr['TEL_AREA'] = self.tel_surface
        hdr['TEL_AREAUNIT'] = 'm2'
        hdr['EXP_TIME'] = self.exp_time
        hdr['EXP_TIMEUNIT'] = 's'

        hdr['ATM_NOISE'] = self.atm_bg_noise
        hdr['ATM_NOISEUNIT'] = 'mag/arcsec2'
        hdr['PHOT_NOISE'] = self.phot_noise
        hdr['RON_NOISE'] = self.readout_noise
        hdr['RON_NOISEUNIT'] = 'e/pix'

        hdr['COMBINE'] = self.combine_method
        hdr['RDI'] = self.RDI
        hdr['ADI'] = self.ADI

        #
        primary_hdu = fits.PrimaryHDU(self.synthobs, header=hdr)
        hdul = fits.HDUList([primary_hdu])

        #export
        name = 'synthobs'
        dir = self.Image.dir
        counter = 0
        candidate = os.path.join(dir, name+'.fits')
        while os.path.exists(candidate):
            counter += 1
            candidate = os.path.join(dir, name+'_'+str(counter)+'.fits')
        hdul.writeto(candidate, overwrite=False)
        print('FITS CREATED :', candidate)


#
def apply_atm_bg(I,wavelength,pixelscale,nx,ny,sky_mag):
    '''
    Apply a constant sky background on the image (2D-array) I.

    Image : scene.Image object

    sky_mag : float
        Sky magnitude/arcsec²
    '''

    print('Adding the sky background with a magnitude/arcsec2 of', sky_mag)
    
    #sky intensity
    F_0 = ref_flux(wavelength)
    I_0 = F_0 / (nx*ny*pixelscale**2)           #intensity in W.m-2.arcsec-2

    I_sky_arcsec2 = I_0*10**(-sky_mag/2.5)             #in W.m-2.arcsec-2
    I_sky = I_sky_arcsec2 * pixelscale**2         #in W.m-2.px-1, unit of I

    #creation of the array and addition to intensity
    I_sky = I_sky*np.ones((nx,ny))
    
    return (I + I_sky)



def apply_photon_noise(I_list,freq,telescope_surface,exp_time,nb_frames,nx,ny):
    '''
    Applies photon noise on an image I in W.m⁻².px⁻¹.

    #args
    Image : scene.Image object
            W.m⁻².px⁻¹

    telescope_surface : float
        The telescope surface in m²
        
    exp_time : float
        Exposure time in s

    '''

    print('Applying the photon noise...')
    
    for frame in range(nb_frames):
        N_ph = Wm2_to_ph(I_list[frame,:,:], nu=freq, surface=telescope_surface, exp_time=exp_time)

        #generation of the noise
        N_noise = np.zeros((nx,ny))
        for i in range(nx):
            for j in range(ny):
                coeff = random.random()
                sign = random.random()
                if sign >= 0.5 :
                    N_noise[i][j] = np.sqrt(N_ph[i][j])*coeff
                else :
                    N_noise[i][j] = -np.sqrt(N_ph[i][j])*coeff
                    if np.abs(N_noise[i][j]) > N_ph[i][j] :   #if noise is over the signal, put the signal to 0
                        N_noise[i][j] = - N_ph[i][j]

        #Addition of the noise
        N_ph = N_ph + N_noise
        I_list[frame,:,:] = ph_to_Wm2(N_ph, nu=freq, surface=telescope_surface, exp_time=exp_time)
    
    return I_list



def apply_readout_noise(I_list,freq,telescope_surface,exp_time,RON,nb_frames,nx,ny):
    '''
    Applies readout noise on an image I in W.m⁻².px⁻¹.

    #args
    Image : scene.Image object
            W.m⁻².px⁻¹

    telescope_surface : float
        The telescope surface in m²
        
    exp_time : float
        Exposure time in s
        
    RON : int
        Readout noise in e/pix

    '''
    
    print('Applying the read-out noise...')
    #
    for frame in range(nb_frames):
        N_ph = Wm2_to_ph(I_list[frame,:,:], nu=freq, surface=telescope_surface, exp_time=exp_time)

        #generation of the noise
        N_noise = np.zeros((nx,ny))
        for i in range(nx):
            for j in range(ny):
                coeff = random.random()
                sign = random.random()
                if sign >= 0.5 :
                    N_noise[i][j] = RON*coeff
                else :
                    N_noise[i][j] = -RON*coeff
                    if np.abs(N_noise[i][j]) > N_ph[i][j] :   #if noise is over the signal, put the signal to 0
                        N_noise[i][j] = - N_ph[i][j]

        #Addition of the noise
        N_ph = N_ph + N_noise
        I_list[frame,:,:] = ph_to_Wm2(N_ph, nu=freq, surface=telescope_surface, exp_time=exp_time)
    
    return I_list
    
        
        

def convolution(I,PSF_on,PSF_off,freq,telescope_surface,exp_time,RON,nb_frames,nx,ny,coro,pixelscale,ADI,angle_list):
    '''
    Take a 2d-single array, mask the central parts of it, convolve it with a PSF sequence and finally add the coronagraph diffraction. Following PSF.need_replication, the on-axis and off-axis PSF have the same nb_frames.

    #args
    Image : scene.Image object
            W.m⁻².px⁻¹
            
    PSF_onaxis : PSF.PSF object
        The on-axis PSF.

    PSF_offaxis : PSF.PSF object
        The off-axis PSF.
        
    pixelscale : float
        Image pixel scale in arcsec
        
    conv_kernel : 2d-array or list
        The PSF sequence. Can be passed as a 2d array if there is only one psf
    coro_diff : 2d-array or list
        Coronagraph diffraction (on-axis PSF) sequence. Can be passed as a 2d array if there is only one psf. If None, must be passed as an empty list.
        
    coronagraph : float or None
        If not None, apply a central opaque mask of radius to be given in mas.
        
    adi_prep : bool, optional
        Prepare the cube for ADI by rotating the MCFOST image before convolving.
        
    angle_list : list or 1d-array
        [only used if adi_prep is True]
        Parallactic angle list for the image to be rotated

    #return
    The convolved sequence
    '''

    #
    print('Convolution...')
    if ADI:
        print('... and rotation')
    
    I_star = copy.deepcopy(np.max(I))  #assuming the intensity of the central star is the one of the brightest pixel
 
    #Central mask
    if coro > 0:
        posx = np.linspace(-nx/2, nx/2, nx)
        posy = np.linspace(-ny/2, ny/2, ny)
        meshx, meshy = np.meshgrid(posx, posy)
        radius_pixel = np.sqrt(meshx ** 2 + meshy ** 2)
        radius_as = radius_pixel * pixelscale
        I[radius_as < coro] = 0.0
    
    #Convolving image and adding diffraction for each psf of the list
    I_list = np.zeros((nb_frames,nx,ny))
    I_mem = copy.deepcopy(I)   #keeping MCFOST image in memory in case of rotation
    for frame in range(nb_frames):

        #rotation for further ADI
        if ADI :
            I = frame_rotate(I_mem, angle_list[frame], imlib='opencv')

        #convolution with the psf
        I_offaxis = np.zeros((nx,ny))
        if PSF_off is not None:
            I_offaxis = convolve_fft(I, PSF_off[frame,:,:])

        #construction of the coronagraphic image
        I_onaxis = np.zeros((nx,ny))
        if PSF_on is not None:
            I_onaxis = I_star * PSF_on[frame,:,:]
        
        #Final image
        I_temp = I_onaxis + I_offaxis

        #put the result in the science image cube
        I_list[frame,:,:] = copy.deepcopy(I_temp)

    return I_list


def recombine(I_list,mode,nb_frames):
    '''
    Recombine an image sequence into a single image.

    #args
    I_list : list-like of 2d-arrays
        The image sequence
    mode : str, optional
        The recombination method in 'median' or 'average' or 'mean'

    #return
    The array that results the recombination
    '''
    
    if nb_frames == 1 :    #if the cube is a single image
        print('The cube is a single image, no need of collapsing.')
    else:
        print('Collapsing the cube. Mode '+mode.lower())
        
    if mode.lower() == 'mean' or mode.lower() == 'average' :
        I_synth = np.mean(I_list, axis=0)
        return I_synth

    elif mode.lower() == 'median' :
        print('Collapsing the cube. Mode '+mode.lower())
        I_synth = np.median(I_list, axis=0)
        return I_synth

    else:
        raise AttributeError('please select a valid combination method : None, mean/average, median')


def RDI(I_list,I_list_RDI,nb_frames,nx,ny):
    '''
    Subtract a reference sequence from an image sequence and returns the result. Both have to be already convolved. (RDI ready data)

    #args
    I_list : nd-array
        The image sequence.
    I_list_RDI : nd-array
        The reference sequence.

    #return
    The RDI sequence
    '''

    #
    print('Reference Differential Imaging')

    #checks
    #if len(I_list) != len(I_list_RDI) :
    #    raise ValueError('Sequences must have the same size ! Here image sequence has '+str(len(I_list))+' exposures while reference sequence has '+str(len(I_list_RDI)))
    if I_list[0].shape != I_list_RDI[0].shape :
        raise ValueError('2d-arrays of the sequences must have the same size ! Here images are '+str(I_list.shape)+' while references are '+str(I_list_RDI.shape))
    
    #reference subtraction
    I_list_sub = np.zeros((nb_frames,nx,ny))
    I_ref_med = np.median(I_list_RDI, axis=0)
    alpha = 1.00
    for frame in range(nb_frames):
        I_list_sub[frame,:,:] = I_list[frame]-alpha*I_ref_med

    return(I_list_sub)



def ADI(I_list, angle_list,nb_frames,nx,ny):
    '''
    I_list :
    
    angle_list :
    '''
    print('Angular Differential Imaging')
    
    #Reference image built from median
    I_subm = np.zeros((nb_frames,nx,ny))
    I_subm_derot = np.zeros((nb_frames,nx,ny))
    I_ref = np.median(I_list, axis=0)
    for frame in range(nb_frames):
        #refine subtraction - TBD
        alpha = 1.00
                
        #Median subtracted cube
        I_subm[frame,:,:] = I_list[frame] - alpha*I_ref
            
        #Derotation
        I_subm_derot[frame,:,:] = frame_rotate(I_subm[frame,:,:], -angle_list[frame])

    return I_subm_derot
