# src/PSF.py
from astropy.io import fits
from scipy.ndimage import zoom
import numpy as np
import copy

class PSF:

    def __init__(self, file=None, onaxis=None, pixelscale=None, FoV=None, nb_frames=None, norm_cst=None, scene_nx=None, scene_pixelscale=None, RDI=False):
        '''file : str
            Path to the  FITS file.

        onxis : boolean, optional
            On-axis or off-axis ?
            
        FoV : float, optional
            Field of view of the instrument.

        nb_frames : int, optional
            Number of frames from the sequence used to make the synthetic observation.
            
        norm_cst : float, optional
            Normalisation constant for the PSF.
            
        '''
    
        #
        self.file = file
        self.onaxis = onaxis
        
        #
        self.pixelscale  = pixelscale
        self.FoV  = FoV
        self.nb_frames = nb_frames
        self.norm_cst = norm_cst
        
        #
        self.scene_nx = scene_nx   #normally the image has been squared
        self.scene_pixelscale = scene_pixelscale
        
        #
        self.RDI = RDI
        
        # Read model results
        self._read()
        
        if file is not None:
            # Print infos
            self.printinfos()
            
            # Check if 3D
            self.check_3D()

            # Select
            self.select()
            
            # Resize - normally the scene image is already square, so no need to check if PSF is square
            self.resize()
            
            # Crop
            self.crop()
            
            # Normalise
            self.normalise()

    def _read(self):
        try:
            if self.file is not None:
                with fits.open(self.file) as hdul:
                    self.image  = hdul[0].data
                    self.header = hdul[0].header
                    self.nx     = hdul[0].header['NAXIS1']
                    self.ny     = hdul[0].header['NAXIS2']
            else:
                self.image  = None
                self.header = None
                self.nx     = None
                self.ny     = None
        except OSError:
            print(f"Cannot open FITS file: {self.file}")
        except Exception as e:
            print(f"Error reading {self.file}: {e}")


    def printinfos(self):
        print('PSF dimensions : ('+str(self.nx)+', '+str(self.ny)+')')
        
    def check_3D(self):
        '''Check that the PSF sequence has 3D (time, x, y), and make it 3D otherwise.'''
        
        #on-axis
        nb_dim = len(self.image.shape)
        if nb_dim <= 1:
            raise ValueError('Input on-axis PSF has only dimension', nb_dim)
        
        elif nb_dim == 2:
            if self.RDI:
                temp_timedim_size = 2*self.nb_frames  #make a sequence 2 times longer for the code below to work and take the end of the sequence
            else:
                temp_timedim_size = self.nb_frames
            temp = np.zeros((temp_timedim_size,self.image.shape[0],self.image.shape[1]))
            for frame in range(self.nb_frames): #duplicate the PSF as many times as requested
                temp[frame,:,:] = copy.deepcopy(self.image)
            self.image = copy.deepcopy(temp)

        elif nb_dim >= 4:
            temp = np.squeeze(self.image)
            nb_dim_temp = len(temp.shape)
            if nb_dim_temp <= 1:
                raise ValueError('Your PSF was squeezed too much')
            elif nb_dim_temp == 2:
                print('Your PSF was squeezed and extended again')
                temp = np.expand_dims(temp, axis=0)
                self.image = copy.deepcopy(temp)
            elif nb_dim_temp == 3:
                print('Your PSF was squeezed')
                self.image = copy.deepcopy(temp)
            else:
                print('Taking the last 3 dimensions of your input array as the PSF')
                self.image = copy.deepcopy(temp[...,:,:,:])
        
        self.nb_frames_tot = self.image.shape[0]

    def resize(self):
        '''to get the right pixel scale'''
        
        #make the thing square if not
        if self.nx != self.ny:
            print('Making the PSF square')
            size = max(self.nx, self.ny)
            top = (size - self.ny) // 2
            bottom = size - self.ny - top
            left = (size - self.nx) // 2
            right = size - self.nx - left

            square_img = np.zeros((self.nb_frames,size,size))
            for frame in range(self.nb_frames):
                square_img[frame:,:,:] = np.pad(
                    self.image[frame,:,:],
                    pad_width=((top, bottom), (left, right)),
                    mode='constant',
                    constant_values=0
                    )
            
            self.image = copy.deepcopy(square_img)
            if self.nx > self.ny:
                self.ny = self.nx
            else:
                self.nx = self.ny
                
        #
        if self.pixelscale is None and self.FoV is not None:
            print('Getting the pixel scale from the FoV and image size')
            self.pixelscale = self.FoV / self.nx
        
        #
        if self.pixelscale != self.scene_pixelscale:
            print('Mismatch between image and PSF pixel scales')
            print('image pixel scale is:', self.scene_pixelscale, 'arsec/pixel')
            print('PSF pixel scale is:', self.pixelscale, 'arsec/pixel')
            print('Resizing the PSF to make them equal...')
            
            zoom_factor = round(self.pixelscale / self.scene_pixelscale)
            
            new_PSF = np.zeros((self.nb_frames,zoom_factor*self.nx,zoom_factor*self.ny))
            for frame in range(self.nb_frames):
                new_PSF[frame,:,:] = zoom(self.image[frame,:,:], zoom_factor, order=1)
            
            new_nx = new_PSF.shape[-1]
            print('New PSF dimensions : ('+str(new_nx)+', '+str(new_nx)+')')
            if self.onaxis and new_nx < self.scene_nx:
                raise ValueError('On-axis PSF is too small compated to the image.')
            
            self.image = copy.deepcopy(new_PSF)
            self.nx = new_nx
            self.ny = new_nx  #PSF is square
            self.pixelscale = self.scene_pixelscale

            
    def crop(self):
        '''if PSF is too large'''
        
        crop_factor = 0.
        
        if self.onaxis and self.nx > self.scene_nx:
            crop_factor = 1.
            print('Cropping the on-axis PSF to', 1/crop_factor, 'of the image size')
            
        elif self.nx > self.scene_nx:
            crop_factor = 1.
            print('Cropping the off-axis PSF to', 1/crop_factor, 'of the image size')
            
        if crop_factor != 0.:
            cx_psf, cy_psf = self.nx//2 + self.nx%2, self.ny//2 + self.ny%2
            print('PSF center is found in coordinate :', cx_psf, cy_psf)
            
            #adjust crop_factor
            offcrop_x, offcrop_y = self.scene_nx//2, self.scene_nx//2  #image is square
            rest_x, rest_y = self.scene_nx%2, self.scene_nx%2
            
            #Centering
            x_min = int(round( (cx_psf -  offcrop_x)        /crop_factor ))
            x_max = int(round( (cx_psf +  offcrop_x+rest_x) /crop_factor ))
            y_min = int(round( (cy_psf -  offcrop_y)        /crop_factor ))
            y_max = int(round( (cy_psf +  offcrop_y+rest_y) /crop_factor ))
        
            if x_min < 0 or y_min < 0 :
                print('WARNING : Cropping too harsh, may not work')
            
            #Slicing
            temp = np.zeros((self.nb_frames, x_max-x_min, y_max-y_min))
            for frame in range(self.nb_frames):
                temp[frame] = self.image[frame, x_min:x_max, y_min:y_max]
            self.image = copy.deepcopy(temp)
        
            self.nx = self.image.shape[-2]
            self.ny = self.image.shape[-1]
            print('New PSF dimensions after cropping : ('+str(self.nx)+', '+str(self.ny)+')')
        else:
            print('No need to crop the PSF')


    def select(self):
        '''according to the number of frames'''
        
        if self.nb_frames_tot < self.nb_frames:
            raise ValueError('You are trying to use more frames than there are in the sequence:', self.nb_frames_tot, 'vs', self.nb_frames )
            
        if self.RDI: #RDI, select frames at the end of the sequence
            if self.nb_frames == 1:
               nb_frames_RDI = 0  #take the first image because there is only one
            else:
                nb_frames_RDI = self.nb_frames_tot - self.nb_frames
                if nb_frames_RDI < self.nb_frames:
                    raise ValueError('Sequence not long enough for RDI, try to use less frames.')
            temp = self.image[nb_frames_RDI:,:,:]   #check that it has the same length science psf
        else: #science, select at the beginning of the sequence
            temp = self.image[:self.nb_frames,:,:]
        self.image = copy.deepcopy(temp)


    def normalise(self):
        '''Normalise to the total flux'''
        
        if self.norm_cst is None:
            self.norm_cst = 1
        print('Normalizing the PSF with a factor', self.norm_cst)
        
        for frame in range(self.nb_frames):
            self.image[frame,:,:] = self.norm_cst * self.image[frame,:,:] / np.sum(self.image[frame,:,:])




def need_replication(PSF_onaxis,PSF_offaxis):
    '''If needed, makes two PSFs with the same number of frames by if nb_frames = 1 for one of them
    
    #args
    PSF_onaxis : PSF.PSF object
    
    PSF_offaxis : PSF.PSF object
    
    #returns
    PSF_onaxis,PSF_offaxis : PSF.PSF objects with the same number of frames.
    '''

    if (PSF_onaxis is None) or (PSF_offaxis is None) :
        print('No need of PSF replication, one of the PSF is None')
        return PSF_onaxis, PSF_offaxis
        
    nb1, nb2 = PSF_onaxis.nb_frames, PSF_offaxis.nb_frames
    if nb1 == nb2:
        print('No need of PSF replication.')
        return PSF_onaxis, PSF_offaxis
        
    elif (nb1 != 1 and nb2 != 1):
        raise ValueError('Please provide the same nb_frames for on-axis and off-axis, or set nb_frames=1 for one of the two.')
        
    else:
        if nb1 == 1:
            print('Replicating the on-axis PSF over the sequence.')
            temp = np.zeros(nb2, PSF_onaxis.nx, PSF_onaxis.ny)
            for frame in range(nb2):
                temp[frame,:,:] = PSF_onaxis.image[0,:,:]
            PSF_onaxis.image = copy.deepcopy(temp)
            PSF_onaxis.nb_frames = nb2
            
        if nb2 == 1:
            print('Replicating the off-axis PSF over the sequence.')
            temp = np.zeros(nb1, PSF_offaxis.nx, PSF_offaxis.ny)
            for frame in range(nb1):
                temp[frame,:,:] = PSF_offaxis.image[0,:,:]
            PSF_offaxis.image     = copy.deepcopy(temp)
            PSF_offaxis.nb_frames = nb1
        
        return PSF_onaxis, PSF_offaxis
