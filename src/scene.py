# src/scene.py
from pymcfost.parameters import Params, find_parameter_file
from astropy.io import fits
import numpy as np
import scipy.constants as sc

class Image:
    '''Inspired by cpinte/pymcfost/image.py'''

    def __init__(self, file=None, mcfost=False, radmc3d=False, **kwargs):
    '''file : str
            Path to the FITS file.
            
        mcfost : boolean, optional
            if in a directory created by mcfost.

        radmc3d : boolean, optional
            if in a directory created by radmc3d.
            
        wl : float, optional
            Wavelength in [µm] (if not in header).
            
        pixelscale : float, optional
            Pixel scale in [arcsec/pixel] (if not in header).
            
        unit : str, optional
            Unit of the image (if not in header).
    '''
    
    # Correct path if needed
    dir = os.path.normpath(os.path.expanduser(dir))
    self.file = file
    self.dir = file.split('/')[:-1]
    
    #
    self.mcfost  = mcfost
    self.radmc3d = radmc3d
    self.wl = wl
    self.pixelscale = pixelscale
    self.unit = unit
    
    #
    if mcfost and radmc3d:
        raise ValueError('Choose your side between mcfost and radmc3d, it cannot be both...')
        
    # mcfost
    if mcfost:
        para_file = find_parameter_file(dir)
        self.mcfostP = Params(para_file)

    # Read model results
    self._read(**kwargs)
    
    # Print infos
    self.printinfos()


    def _read(self):
        #
        if mcfost:
            try:
                hdu = fits.open(self.dir + "/RT.fits.gz")
                self.image = hdu[0].data
                # Read a few keywords in header
                self.header = hdu[0].header
                self.pixelscale = hdu[0].header['CDELT2'] * 3600.0  # arcsec
                self.unit = hdu[0].header['BUNIT']
                self.wl = hdu[0].header['WAVE']  # micron
                self.freq = sc.c / (self.wl * 1e-6)
                self.cx = hdu[0].header['CRPIX1']
                self.cy = hdu[0].header['CRPIX2']
                self.nx = hdu[0].header['NAXIS1']
                self.ny = hdu[0].header['NAXIS2']
                try:
                    self.star_positions = hdu[1].data
                except:
                    self.star_positions = []
                try:
                    self.star_vr = hdu[2].data
                except:
                    self.star_vr = []
                try:
                    self.star_properties = hdu[3].data
                except:
                    self.star_properties = []
                hdu.close()
            except OSError:
                print('cannot open', self._RT_file)

        #radmc3D
        elif radmc3d:
            pass
        
        #general case
        else:
            try:
                with fits.open(self.filename) as hdul:
                    self.image = hdul[0].data
                    self.header = hdul[0].header

                    # --- Try to read from header; fall back to user-provided values ---
                    self.pixelscale = (
                        self.header.get("CDELT2", None) * 3600.0
                        if "CDELT2" in self.header
                        else self.pixelscale
                    )
                    self.unit = self.header.get("BUNIT", self.unit)
                    self.wl = self.header.get("WAVE", self.wl)

                    # --- Derived quantities ---
                    if self.wl is not None:
                        self.freq = sc.c / (self.wl * 1e-6)  # Hz
                    else:
                        self.freq = None

                    # --- Image shape and center ---
                    self.ny, self.nx = self.image.shape
                    self.cx = self.header.get("CRPIX1", self.nx / 2)
                    self.cy = self.header.get("CRPIX2", self.ny / 2)
                    if cx != int(cx):
                        cx = int(cx) +1
                    if cy != int(cy):
                        cy = int(cy) +1
                        
        except OSError:
            print(f"Cannot open FITS file: {self.filename}")
        except Exception as e:
            print(f"Error reading {self.filename}: {e}")


    def printinfos(self):
    
        print('lambda = ', self.wl, 'µm')
        print('pixel scale = ', self.pixelscale, 'arcsec/pix')
        print('[nx, ny] = ', [self.nx, self.ny], 'pix')
        print('[cx, cy] = ', [self.cx, self.cy], 'pix')
        print('Input image is in', self.unit)
        
    
    def convert(self):
        

    
