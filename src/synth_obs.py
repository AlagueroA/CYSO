# src/synth_obs.py

from .utils import *



def apply_atm_bg(I, sky_mag, wavelength, pixel_scale) :
    '''
    Apply a constant sky background on the image (2D-array) I.

    I : ndarray
        A single 2d-array
    sky_mag : float
        Sky magnitude/arcsec²
    wavelength : float
        Wavelength in µm
    pixel_scale : float
        pixel scale in arcsec
    '''


    #sky intensity
    F_0 = ref_flux(wavelength)
    I_0 = F_0 / (len(I)*len(I[0])*pixel_scale**2)           #intensity in W.m-2.arcsec-2

    I_sky_arcsec2 = I_0*10**(-sky_mag/2.5)             #in W.m-2.arcsec-2
    I_sky = I_sky_arcsec2 * pixel_scale**2         #in W.m-2.px-1, unit of I
    if np.mean(I) > 1e-12 :  #put to Jy if needed
        I_sky = Wm2_to_Jy(I_sky, sc.c/(wavelength*1e-6))
        #print(I_sky)
        #print(np.mean(I))

    #creation of the array and addition to intensity
    I_sky = I_sky*np.ones((len(I),len(I[0])))


    I = I + I_sky
    
    return I



def apply_photon_noise(I, wavelength, telescope_surface, exp_time) :
    '''
    Applies photon noise on an image I in W.m⁻².px⁻¹.

    #args
    I : ndarray
        A single 2d-array in W.m⁻².px⁻¹
    wavelength : float
        Wavelength in µm, needed to convert to a photon number
    telescope_surface : float
        The telescope surface in m²
    exp_time : float
        Exposure time in s

    #return
    The array with photon noise applied

    '''
    
    Jsk = False

    #conversion of I to a photon number
    freq = sc.c / (wavelength*1e-6)
    if np.mean(I) > 1e-12 :  #if perhaps the map is in Jy, convert it to Wm-2
        #print('Jy detected, converting to Wm-2')
        Jsk = True
        I = Jy_to_Wm2(I, freq)
        
    N_ph = Wm2_to_ph(I, nu=freq, surface=telescope_surface, exp_time=exp_time)

    #generation of the noise
    N_noise = np.zeros((len(N_ph),len(N_ph[0])))
    for i in range(len(N_ph)):
        for j in range(len(N_ph[0])):
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
    I = ph_to_Wm2(N_ph, nu=freq, surface=telescope_surface, exp_time=exp_time)
    
    if Jsk : #convert back to Jy if needed
        I = Wm2_to_Jy(I, freq)
        

    return(I)


def apply_readout_noise(I, wavelength, telescope_surface, exp_time, RON=2) :
    '''
    Applies readout noise on an image I in W.m⁻².px⁻¹.

    #args
    I : ndarray
        A single 2d-array in W.m⁻².px⁻¹
    wavelength : float
        Wavelength in µm, needed to convert to a photon number
    telescope_surface : float
        The telescope surface in m²
    exp_time : float
        Exposure time in s
    RON : int
        Readout noise in e/pix

    #return
    The array with readout noise applied

    '''

    #conversion of I to a photon number
    freq = sc.c / (wavelength*1e-6)
    if np.mean(I) > 1e-12 :  #if perhaps the map is in Jy, convert it to Wm-2
        Jsk = True
        I = Jy_to_Wm2(I, freq)
    N_ph = Wm2_to_ph(I, nu=freq, surface=telescope_surface, exp_time=exp_time)

    #generation of the noise
    N_noise = np.zeros((len(N_ph),len(N_ph[0])))
    for i in range(len(N_ph)):
        for j in range(len(N_ph[0])):
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
    I = ph_to_Wm2(N_ph, nu=freq, surface=telescope_surface, exp_time=exp_time)
    
    if Jsk : #convert back to Jy if needed
        I = Wm2_to_Jy(I, freq)

    return(I)


def convolution(I, pixelscale, conv_kernel, coro_diff, coronagraph=None, adi_prep=False, angle_list=None):
    '''
    Take a 2d-single array, mask the central parts of it, convolve it with a PSF sequence and finally add the coronagraph diffraction.

    #args
    I : ndarray
        A single 2d-array
    pixelscale : float
        Image pixel scale in arcsec
    conv_kernel : 2d-array or list
        The PSF sequence. Can be passed as a 2d array if there is only one psf
    coro_diff : 2d-array or list
        Coronagraph diffraction (on-axis PSF) sequence. Can be passed as a 2d array if there is only one psf. If None, must be passed as an empty list.
    coronagraph : float or None
        If not None, apply a central opaque mask of radius to be given in arcsec
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

    nx = len(I)
    ny = len(I[0])

    try :
        if conv_kernel == []:       #if no convolution
            conv_kernel=None
            I_list = [I]
    except ValueError :
        pass

    if conv_kernel is not None :
        if len(conv_kernel.shape) == 2 :           #if psf is a single 2d-array, making it a list
            conv_kernel = [conv_kernel]
        try :
            if coro_diff==[] :                              #if no coronagraph diffraction, initializing an array of 0
                coro_diff = np.zeros((nx, ny))
        except ValueError :
            pass
        if len(coro_diff.shape) == 2 :             #if coro. diffraction is a single 2d-array, making it a list of the same length of the psf list
            coro_diff_temp = copy.deepcopy(coro_diff)
            coro_diff = []
            for i_ in range(len(conv_kernel)):
                coro_diff.append(coro_diff_temp)
        if len(conv_kernel) != len(coro_diff) :
            raise ValueError('PSF list musts be the same size of the coronagraph diffraction list')

        #Lists that will contain the convolved images
        I_list = []

        #Assuming the intensity of the central star is the one of the brightest pixel
        I_star = np.max(I)

        #Central mask
        if coronagraph is not None :

            posx = np.linspace(-nx/2, nx/2, nx)
            posy = np.linspace(-ny/2, ny/2, ny)
            meshx, meshy = np.meshgrid(posx, posy)
            radius_pixel = np.sqrt(meshx ** 2 + meshy ** 2)
            radius_as = radius_pixel * pixelscale

            I[radius_as < coronagraph] = 0.0

        #Convolving image and adding diffraction for each psf of the list
        compt = 0
        I_mem = copy.deepcopy(I)   #keeping MCFOST image in memory inc ase of rotation
        for diff, psf in zip(coro_diff, conv_kernel) :

            #rotation for further ADI
            if adi_prep :
                if compt == 0 :
                    print('... and rotation')
                I = frame_rotate(I_mem, angle_list[compt], imlib='opencv')
                compt += 1

            #convolution with the psf
            I_temp = convolve_fft(I, psf)

            #construction of the coronagraphic image
            I_coro = I_star * diff
            I_temp = I_temp + I_coro

            #put the result in the science image cube
            I_list.append(I_temp)

    return(I_list)



def recombine(I_list, mode='median'):
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

    #
    I_list = np.array(I_list)
    
    if len(I_list.shape) == 2 :    #if the cube is a single image
        I = I_list
        print('The cube is a single image, no need of collapsing.')
        
    elif mode.lower() == 'mean' or mode.lower() == 'average' :
        print('Collapsing the cube. Mode '+mode.lower())
        I = np.zeros((len(I_list[0]), len(I_list[0][0])))
        for ind_x in range(len(I_list[0])):
             for ind_y in range(len(I_list[0][0])):
                 for ind_I in range(len(I_list)):
                    I[ind_x][ind_y] += I_list[ind_I][ind_x][ind_y] / len(I_list)

    elif mode.lower() == 'median' :
        print('Collapsing the cube. Mode '+mode.lower())
        I = np.zeros((len(I_list[0]), len(I_list[0][0])))
        for ind_x in range(len(I_list[0])):
             for ind_y in range(len(I_list[0][0])):
                pix_list_I = []
                for ind_I in range(len(I_list)):
                    pix_list_I.append(I_list[ind_I][ind_x][ind_y])
                I[ind_x][ind_y] = np.median(pix_list_I)

        I = np.median(I_list, axis=0)

    else:
        raise AttributeError('please select a valid combination method : None, mean/average, median')

    return(I)


def RDI(I_list, I_ref_list):
    '''
    Subtract a reference sequence from an image sequence and returns the result. Both have to be already convolved. (RDI ready data)

    #args
    I_list : 2d-array or list
        The image sequence as a list or a 2d-array if it is a single image.
    I_ref_list : 2d-array or list
        The reference sequence. Must be the same shape as I_list.

    #return
    The RDI sequence
    '''

    #
    print('applying RDI')

    #checks
    if not isinstance(I_list, list):
        I_list = [I_list]
    if not isinstance(I_ref_list, list):
        I_ref_list = [I_ref_list]
    if len(I_list) != len(I_ref_list) :
        raise ValueError('Sequences must have the same size ! Here image sequence has '+str(len(I_list))+' exposures while reference sequence has '+str(len(I_ref_list)))
    if I_list[0].shape != I_ref_list[0].shape :
        raise ValueError('2d-arrays of the sequences must have the same size ! Here images are '+str(I_list.shape)+' while references are '+str(I_ref_list.shape))
    
    #
    nx = I_list[0].shape[0]
    
    #reference subtraction
    I_RDI_list = []
    I_ref_med = np.median(np.array(I_ref_list), axis=0)
    alpha = 1.00
    for ind_I in range(len(I_list)):
        #I_RDI_list.append(I_list[ind_I]-I_ref_list[ind_I])
        I_RDI_list.append(I_list[ind_I]-alpha*I_ref_med)
        
#    I_ref_med = np.median(np.array(I_ref_list), axis=0)
#    for ind_I in range(len(I_list)):
#        rest = 1e80
#        alpha = 1
#        nb_neg = nx/2
#        test_alpha = np.linspace(-2, 2, 200)
#        for xfactor in test_alpha :
#            test = np.array(I_list[ind_I])-xfactor*I_ref_med
#            count_neg = np.sum(test<0)
#            rest_test = np.sum( (test)**2 )
#            if rest_test < rest :
#                rest = rest_test
#                alpha = xfactor
#                #nb_neg = count_neg
                
#        I_RDI_list.append(I_list[ind_I]-alpha*I_ref_med)

    return(I_RDI_list)



def ADI(I_list, angle_list) :
    '''
    I_list :
    
    angle_list :
    '''
    
    #Converting into numpy array to be sure
    I_list = np.array(I_list)
    
    #Reference image built from median
    I_subm = []
    I_ref = np.ones(I_list.shape)
    for ind_psf in range(I_list.shape[0]):
        I_ref[ind_psf,:,:] = np.median(I_list, axis=0)
        
        #refine subtraction
        alpha = 0.99
#        rest = 1e80
#        alpha = 1
#        nb_neg = nx/2
#        test_alpha = np.linspace(-2, 2, 200)
#        for xfactor in test_alpha :
#            test = np.array(I_list[ind_psf])-xfactor*np.array(I_ref_list[ind_psf])
#            count_neg = np.sum(test<0)
#            rest_test = np.sum( (test)**2 )
#            if rest_test < rest :
#                rest = rest_test
#                alpha = xfactor
                
        #Median subtracted cube
        I_subm.append(I_list[ind_psf] - alpha*I_ref[ind_psf])
    

            
    #Derotation
    I_subm_derot = np.ones(I_list.shape)
    for ind_psf in range(I_list.shape[0]):
        I_subm_derot[ind_psf] = frame_rotate(I_subm[ind_psf], -angle_list[ind_psf])
    
    #I_tot = np.median(I_subm_derot, axis=0)
    
    return I_subm_derot
    


def SNR(image, mode='map', telescope_diameter=8, image_ref=None, n_patch=10, i=0, iaz=0, i_star=0, i_planet=1, incl=0, planet_pos=None, ignore_neg=False, **kwargs):
    '''
    Signal to Noise Ratio,  output depends on the selected mode.

    image : pymcfost image.Image object

    mode : str, optional
        The recombination method in 'map' or 'planet' or 'diff'

        'planet' : estimate the signal in a fwhm around the planet and the noise in a fwhm-wide annulus at the same separation

        'diff' : estimate the signal in a fwhm around the planet on I and the noise at the same location in I_ref
            I_ref : 2d-array
                Intensity map from the same simulation as I but without the planet. Should not be None if diff.

        'patch' : estimate the signal in a fwhm around the planet and the noise as the rms of 2fwhm-wide patches at the same separation. SNR calculated as described in Mawet+2014.
            n_patch : int, optional
                Total number of patches on which the noise is estimated. Min is 2.

        'map' : plot an SNR map following vip_hci,kw  arguments can be put as kwargs
            array : numpy ndarray
                Input frame (2d array). Is the intensity in image.
            fwhm : float
                Size in pixels of the FWHM. Calculated from image.
            approximated : bool, optional
                If True, an approximated S/N map is generated.
            plot : bool, optional
                If True plots the S/N map. True by default.
            known_sources : None, tuple or tuple of tuples, optional
                To take into account existing sources. It should be a tuple of float/int
                or a tuple of tuples (of float/int) with the coordinate(s) of the known
                sources.
            nproc : int or None
                Number of processes for parallel computing.
            array2 : numpy ndarray, optional
                Additional image (e.g. processed image with negative derotation angles)
                enabling to have more noise samples. Should have the
                same dimensions as array. Taken from image_ref.
            use2alone: bool, optional
                Whether to use array2 alone to estimate the noise (might be useful to
                estimate the snr of extended disk features).
            verbose: bool, optional
                Whether to print timing or not.
            **kwargs : dictionary, optional
                Arguments to be passed to ``plot_frames`` to customize the plot (and to
                save it to disk).

    telescope_diameter : float, optional
        Diameter of the telescope primary mirror in m. Default is VLT.

    image_ref : pymcfost Image.image object, optional
        Reference image for the noise to evaluate the noise. Needed for mode='diff' and eventually for mode='map'.

    i : int, optional
        Inclination index of the mcfost image. Default is 0.

    i_az : int, optional
        Azimuth index of the mcfost image. Default is 0.

    i_star : int, optional
        Star index in the mcfost star list. Default is 0.

    i_planet : int, optional
        Planet index in the mcfost star list. Default is 1.
        
    incl : float, optional
        Inclination of the system in rad

    planet_pos : tuple, optional
        (x,y) coordinates of the planet in pixels
        
    **kwargs : for hci_vip.snrmap

    '''

    #Intensity
    try :
        I = image.image[0, iaz, i, :, :]
        I_ref = None
        if image_ref is not None :
            I_ref = image_ref.image[0, iaz, i, :, :]
    except IndexError :
        I = image.image[iaz, i, :, :]
        I_ref = None
        if image_ref is not None :
            I_ref = image_ref.image[iaz, i, :, :]
        

    #fwhm
    fwhm_pix = (((image.wl*1e-6)/telescope_diameter) * 180*3600/np.pi) / image.pixelscale
    fwhm_pix_ = 6

    #planet position
    if planet_pos is not None :
        x_planet, y_planet = planet_pos[0], planet_pos[1]
    else :
        x_planet = int(-image.star_positions[0][0][0][i_planet]/image.pixelscale)    #in pixels
        y_planet = int(image.star_positions[1][0][0][i_planet]/image.pixelscale)
        
    r_planet = np.sqrt( (x_planet)**2 + (y_planet)**2 )   #true radius taking into account inclination
    phi_planet = np.arctan2(y_planet, x_planet)
    
    #separation grids
    halfsize = np.asarray(image.image.shape[-2:]) / 2
    posx = np.linspace(-halfsize[0], halfsize[0], image.nx)
    posy = np.linspace(-halfsize[1], halfsize[1], image.ny)
    meshx_, meshy_ = np.meshgrid(posx, posy)
    meshx, meshy = np.meshgrid(posx*np.cos(incl), posy)
    r_grid = np.sqrt((meshx) ** 2 + (meshy) ** 2)
    r_grid_planet = np.sqrt((meshx_-x_planet) ** 2 + (meshy_-y_planet) ** 2) #np.sqrt((meshx_-x_planet-3) ** 2 + (meshy_-y_planet-3) ** 2)   #modify this to fit planet position on the inclined grid
    r_grid_corner = np.sqrt((meshx+x_planet) ** 2 + (meshy+y_planet) ** 2)
    phi_grid = np.arctan2(meshy, meshx) #origin on the right
    
    #masks
    mask_planet = r_grid_planet <= 2*fwhm_pix     #mask around planet
    mask_annulus = np.logical_and(r_planet*np.cos(incl) - 1*fwhm_pix <= r_grid, r_grid <= r_planet*np.cos(incl) + 1*fwhm_pix)      #inclined annulus
    mask_noise = np.logical_and(mask_annulus, np.logical_not(mask_planet))

    if ignore_neg :
        mask_planet = np.logical_and(mask_planet, I>0)
        mask_annulus = np.logical_and(mask_annulus, I>0)
        mask_noise = np.logical_and(mask_noise, I>0)
        
        
    #patches - around the planet
    incr_angle = 2*np.pi/(n_patch)     #add   (n_patch+1) to remove the patch on the right
    off_x, off_y = 5*fwhm_pix_, 0*fwhm_pix_
    x_patch, y_patch = (x_planet + off_x) * np.cos(incl), y_planet + off_y
    masks_patch = []
    for ind in range(n_patch):
#        x_patch_old, y_patch_old = x_patch, y_patch
#        x_patch = np.cos(incr_angle)*x_patch_old - np.sin(incr_angle)*y_patch_old
#        y_patch = np.sin(incr_angle)*x_patch_old + np.cos(incr_angle)*y_patch_old
        off_x_old, off_y_old = off_x, off_y
        off_x = np.cos(incr_angle)*off_x_old - np.sin(incr_angle)*off_y_old
        off_y = np.sin(incr_angle)*off_x_old + np.cos(incr_angle)*off_y_old
        x_patch, y_patch = (x_planet + off_x) * np.cos(incl), y_planet + off_y
        x_patch_ = x_patch / np.cos(incl)
        
        r_grid_patch = np.sqrt( ((meshx_-x_patch_))**2 + ((meshy_-y_patch))**2 )   #round patches along an ellipse in a res elmt
        if ignore_neg :
            mask_patch = np.logical_and(r_grid_patch <= 1.22*fwhm_pix_, I>=0)
        else :
            mask_patch = r_grid_patch <= 1.22*fwhm_pix_
        masks_patch.append(mask_patch)

    #creation of the final noise mask - removing the first mask of the list, ie the planet
    mask_patch = masks_patch[0]
    for ind, mask in enumerate(masks_patch):
       # if ind != 0 :
        mask_patch = np.logical_or(mask_patch, mask)
    #plt.imshow(np.logical_or(mask_planet, mask_patch), origin='lower')
    #plt.imshow(mask_planet, origin='lower')
    #plt.imshow(mask_patch, origin='lower')
    #plt.imshow(I, alpha=0.5, cmap='twilight',origin='lower')
    #plt.show()
    #exit()

    if mode=='planet':

        print('Calculating the SNR of the planet...')
        I_signal = np.where(mask_planet, I, np.nan)
        I_noise = np.where(mask_noise, I, np.nan)
        

        #S = np.sum(I_signal)    #signal in a resolution element around the planet
        S = np.nanmax(I_signal)    #signal in as peak as the planet position
        #N = np.sqrt(np.median(I_noise**2)) #* len(I_signal)  #noise = rms of the annulus * nb of pixels inside a resolution element
        N = np.nanmedian(I_noise) #noise = median of the annulus * nb of pixels inside a resolution element
        
        if N < 0 and S > 0: #handle negative noise
            SNR = ( S + 2*np.abs(N) ) / np.abs(N)
        else :
            SNR = S/N

        print(S, N)

        delta_S = 0
        delta_N = 0
        delta_SNR = 0#np.sqrt( (delta_S/N)**2 + (delta_N * S/N**2)**2 )

        #print('SNR =', SNR, ' ; delta_SNR = ', delta_SNR)
        #return(SNR, delta_SNR)
        return SNR
        

    if mode=='patch':
        print('Calculating SNR of the planet from patches...')
        F_signal = np.nanmedian(np.where(mask_planet,I,np.nan)) #np.sum(I*np.where(I,mask_planet,0))
        F_patch = np.nanmedian(np.where(mask_patch,I,np.nan))
        std_patch = np.nanstd(np.where(mask_patch,I,np.nan))
#        F_patch = []
#        for patch in masks_patch :
#            F_patch.append(np.nanmean(np.where(patch,I,0)))
#        F_patch = np.array(F_patch)

        #S = F_signal - np.mean(F_patch)
        #N = np.std(F_patch)*np.sqrt(1+1/n_patch)
        #Mawet+2014
        S = F_signal                      #x1
        N = F_patch #np.nanmean(F_patch)              #x2
        delta_N = std_patch#np.nanstd(F_patch)         #s2
        

        SNR = (S - N) / (delta_N * np.sqrt(1+1/n_patch))
        #SNR = S / N
        print('SNR =', SNR)
        return(SNR)

    if mode=='corner':
        print('Calculating SNR of the planet...')
        print('Evaluating noise in the corner of the image...')
        
        #TBD
    
    if mode=='annulus': #calculate the signal in an annulus and compare to the signal at the position of the planet
        print('Annulus mode detected, trying to find the peak created by the planet...')
        
        plt.imshow(np.where(mask_annulus, I, 0), origin='lower')
        plt.show()
        
        az_bins = np.linspace(-np.pi, np.pi, 360)   #bins of 1° width
        I_az = np.zeros(az_bins.shape[0]-1)
        for j in range(az_bins.shape[0]-1):
            mask_az = np.logical_and(az_bins[j]<phi_grid, phi_grid<az_bins[j+1])
            mask_loc = np.logical_and(mask_annulus, mask_az)
            I_az[j] = np.sum( np.where(mask_loc, I, 0) )
        
        I_az = savgol_filter(I_az, window_length=20, polyorder=2)
        plt.plot(az_bins[:-1], I_az)
        plt.axvline(x=az_bins[:-1][np.searchsorted(az_bins[:-1], phi_planet)] )
        

        def gaussian_off(x, A, x0, sigma, coeff, b):
            return A * np.exp( - (x-x0)**2 / (2*sigma**2) ) + coeff*x + b
        
        fwhm_pix_az = 2* np.tan( fwhm_pix / r_planet )
        sigma_pix_az = fwhm_pix_az / ( 2.0 * np.sqrt(2.0 * np.log(2)) )
        def model_gaussian_off(x, A, sigma, coeff, b):
            #return gaussian_off(x=x, A=A, x0=phi_planet, sigma=fwhm_pix_az, coeff=coeff, b=b)
            return gaussian_off(x=x, A=A, x0=phi_planet, sigma=sigma, coeff=coeff, b=b)

        ind_planet = np.searchsorted(az_bins[:-1], phi_planet)
        az_bins_fit = az_bins[:-1][ind_planet-50:ind_planet+50]   #100° around the planet
        I_az_fit = I_az[ind_planet-50:ind_planet+50]
        
        popt, pcov = curve_fit(model_gaussian_off, az_bins_fit, I_az_fit)
        
        plt.plot(az_bins[:-1], model_gaussian_off(az_bins[:-1], *popt))
        plt.show()
        
        #check if the fit is real or not : if sigma is larger than the PSF, the fit is wrong
        if popt[1] > 10*sigma_pix_az or popt[0] < 0 :
            print('Planet not detected')
            SNR = 0
            return SNR
        
        SNR = popt[0] / (popt[2]*phi_planet+popt[3])    # = A / (a*x0+b)
        print('SNR of the planet :', SNR)
        
        return SNR
        
    if mode=='map':

        print('Using vip_hci snrmap : ')

        snrmap(I, fwhm=fwhm_pix, plot=True, array2=I_ref, **kwargs)



def SNR_disk(image, i=0, iaz=0, nb_rings=100, r_min=1, r_max=200, incl=0, r_noise=None, return_r=False, phot_noise=False):
    '''
    Do if disk is inclined.
    '''
    
    #Intensity
    I = image.image[0, iaz, i, :, :]
    
    #Grids
    halfsize = np.asarray(image.image.shape[-2:]) / 2
    posx = np.linspace(-halfsize[0], halfsize[0], image.nx) * image.pixelscale
    posy = np.linspace(-halfsize[1], halfsize[1], image.ny) * image.pixelscale
    meshx, meshy = np.meshgrid(posx, posy)
    r_grid = np.sqrt( (meshx*np.cos(incl)) ** 2 + (meshy) ** 2)
    
    #Signal
    rings = np.linspace(r_min, r_max, nb_rings)
    S = np.zeros(nb_rings-1)
    if phot_noise :
        N = np.zeros(nb_rings-1)
    for ind in range(nb_rings-1):
        mask = np.logical_and(rings[ind]<r_grid, r_grid<rings[ind+1])
        S[ind] = np.median( I[mask] )
        N[ind] = np.sqrt( S[ind] )
        
        
    #Noise
    if not phot_noise :
        if r_noise is None :
            r_noise = 2*r_max
        r_min_noise, r_max_noise = 0.9*r_noise, 1.1*r_noise
        mask_noise = np.logical_and(r_min_noise<r_grid, r_grid<r_max_noise)
        N = np.sqrt( np.mean( I[mask_noise]**2 ) )
    
#    plt.figure()
#    plt.imshow(I, origin='lower')
#    plt.show()

    SNR = S/N

    
    if return_r :
        return SNR, rings
    else :
        return SNR
    


def r_profile(image, i=0, iaz=0, nb_rings=100, r_min=1, r_max=200, incl=0, top_half=False, r_noise=None, return_r=False, ret_rms=False):
    '''
    Do if disk is inclined.
    '''
    
    #Intensity
    try :
        I = image.image[0, iaz, i, :, :]
    except IndexError :
        I = image.image[iaz, i, :, :]
    
    #Grids
    halfsize = np.asarray(image.image.shape[-2:]) / 2
    posx = np.linspace(-halfsize[0], halfsize[0], image.nx) * image.pixelscale
    posy = np.linspace(-halfsize[1], halfsize[1], image.ny) * image.pixelscale
    meshx, meshy = np.meshgrid(posx, posy)
    r_grid = np.sqrt( (meshx*np.cos(incl)) ** 2 + (meshy) ** 2)
    if top_half :
        mask_half = np.zeros_like(I, dtype=bool)
        mask_half[int(halfsize[1]):, :] = True
    
    #Signal
    rings = np.linspace(r_min, r_max, nb_rings)
    S = np.zeros(nb_rings-1)
    N = np.zeros(nb_rings-1)
    for ind in range(nb_rings-1):
        mask = np.logical_and(rings[ind]<r_grid, r_grid<rings[ind+1])
        if top_half :
            mask = np.logical_and(mask, mask_half)
        S[ind] = np.nanmedian( I[mask] )#np.sum( I[mask] )#np.mean(I, where=mask)   #np.sum( I[mask] )#np.median( I[mask] )
        if ret_rms :
            S[ind] = np.sqrt( np.nanmean( I[mask]**2 ) )
        N[ind] = np.sqrt(np.sum( I[mask] ) )#np.median( I[mask] )
    
    #plt.imshow(mask, origin='lower')
    #plt.show()
    
    if return_r :
        return S, N, rings
    else :
        return S

