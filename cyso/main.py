# src/main.py
#import .CYSO_input as ui  #to read from the active directory, user needs to copy it
from pathlib import Path
import runpy

from .scene import *
from .PSF import *
from .synth_obs import *

## READ INPUT
#path_image           = ui.path_image
#mcfost               = ui.mcfost
#radmc3D              = ui.radmc3D
#wavelength           = ui.wavelength
#pixelscale           = ui.pixelscale
#unit                 = ui.unit
#distance             = ui.distance
#path_pa              = ui.path_pa
#lat                  = ui.lat
#deg                  = ui.deg
#coro                 = ui.coro
#tel_diam             = ui.tel_diam
#tel_surface          = ui.tel_surface
#phot_noise           = ui.phot_noise
#readout_noise        = ui.readout_noise
#atm_bg_noise         = ui.atm_bg_noise
#path_PSF_onaxis      = ui.path_PSF_onaxis
#path_PSF_offaxis     = ui.path_PSF_offaxis
#nb_frames            = ui.nb_frames
#exp_time             = ui.exp_time
#combine_method       = ui.combine_method
#pixelscale_PSF       = ui.pixelscale_PSF
#FoV                  = ui.FoV
#norm_cst             = ui.norm_cst
#RDI                  = ui.RDI
#path_image_RDI       = ui.path_image_RDI
#path_PSF_onaxis_RDI  = ui.path_PSF_onaxis_RDI
#path_PSF_offaxis_RDI = ui.path_PSF_offaxis_RDI
#ADI                  = ui.ADI
#export_fits          = ui.export_fits


def run_main():
    ## READ INPUT
    if not Path('CYSO_input.py').exists():
        print('CYSO_input.py not found in current directory. Run cyso --make_setup first.')
        return

    # import user inputs
    ui = runpy.run_path('CYSO_input.py')
    print('Running cyso with CYSO_input.py contents:')
    #print(user_data)  #find a better way to plot the arguments loaded
    #exit()
    
    path_image           = ui['path_image']
    mcfost               = ui['mcfost']
    radmc3d              = ui['radmc3d']
    wavelength           = ui['wavelength']
    pixelscale           = ui['pixelscale']
    unit                 = ui['unit']
    distance             = ui['distance']
    path_pa              = ui['path_pa']
    lat                  = ui['lat']
    dec                  = ui['dec']
    coro                 = ui['coro']
    tel_diam             = ui['tel_diam']
    tel_surface          = ui['tel_surface']
    phot_noise           = ui['phot_noise']
    readout_noise        = ui['readout_noise']
    atm_bg_noise         = ui['atm_bg_noise']
    path_PSF_onaxis      = ui['path_PSF_onaxis']
    path_PSF_offaxis     = ui['path_PSF_offaxis']
    nb_frames            = ui['nb_frames']
    exp_time             = ui['exp_time']
    combine_method       = ui['combine_method']
    pixelscale_PSF       = ui['pixelscale_PSF']
    FoV                  = ui['FoV']
    norm_cst             = ui['norm_cst']
    RDI                  = ui['RDI']
    path_image_RDI       = ui['path_image_RDI']
    path_PSF_onaxis_RDI  = ui['path_PSF_onaxis_RDI']
    path_PSF_offaxis_RDI = ui['path_PSF_offaxis_RDI']
    ADI                  = ui['ADI']
    export_fits          = ui['export_fits']
    
    
    #
    if RDI:
        print('\nYou put RDI=True, verbose will appear in double.')
    

    ## IMAGE
    print('\n------------ IMAGE ------------')
    print('Loading the image...')

    image = Image(file=path_image, mcfost=mcfost, radmc3d=radmc3d, wl=wavelength, pixelscale=pixelscale, unit=unit)
    
    image_RDI = None
    if RDI:
        image_RDI = Image(file=path_image_RDI, mcfost=mcfost, radmc3d=radmc3d, wl=wavelength, pixelscale=pixelscale, unit=unit)

    ## PSF
    print('\n------------ PSF ------------')
    print('Loading the PSF and preparing it like a chef...')

    PSF_onaxis = PSF(file=path_PSF_onaxis,
                     onaxis=True,
                     pixelscale=pixelscale_PSF,
                     FoV=FoV,
                     nb_frames=nb_frames,
                     norm_cst=norm_cst,
                     scene_nx=image.nx,
                     scene_pixelscale=image.pixelscale,
                     RDI=False)

    PSF_offaxis = PSF(file=path_PSF_offaxis,
                     onaxis=False,
                     pixelscale=pixelscale_PSF,
                     FoV=FoV,
                     nb_frames=nb_frames,
                     norm_cst=norm_cst,
                     scene_nx=image.nx,
                     scene_pixelscale=image.pixelscale,
                     RDI=False)

    PSF_onaxis, PSF_offaxis = need_replication(PSF_onaxis, PSF_offaxis)

    PSF_onaxisRDI, PSF_offaxisRDI = None, None
    if RDI:
        PSF_onaxisRDI = PSF(file=path_PSF_onaxis,
                     onaxis=True,
                     pixelscale=pixelscale_PSF,
                     FoV=FoV,
                     nb_frames=nb_frames,
                     norm_cst=norm_cst,
                     scene_nx=image.nx,
                     scene_pixelscale=image.pixelscale,
                     RDI=True)

        PSF_offaxisRDI = PSF(file=path_PSF_offaxis,
                     onaxis=False,
                     pixelscale=pixelscale_PSF,
                     FoV=FoV,
                     nb_frames=nb_frames,
                     norm_cst=norm_cst,
                     scene_nx=image.nx,
                     scene_pixelscale=image.pixelscale,
                     RDI=True)

        PSF_onaxisRDI, PSF_offaxisRDI = need_replication(PSF_onaxisRDI, PSF_offaxisRDI)
        
        

    ## SYNTHETIC OBSERVATION
    print('\n------------ SYNTHETIC OBS. ------------')
    print('Now cooking the synthetic observation...')

    synthobs = SynthObs(Image=image, PSF_onaxis=PSF_onaxis, PSF_offaxis=PSF_offaxis,
                        coro=coro,
                        tel_diam=tel_diam,
                        tel_surface=tel_surface,
                        phot_noise=phot_noise,
                        exp_time=exp_time,
                        readout_noise=readout_noise,
                        atm_bg_noise=atm_bg_noise,
                        combine_method=combine_method,
                        RDI=RDI,
                        Image_RDI=image_RDI,
                        PSF_onaxisRDI=PSF_onaxisRDI,
                        PSF_offaxisRDI=PSF_offaxisRDI,
                        ADI=ADI,
                        path_pa=path_pa,
                        lat=lat,
                        dec=dec,
                        export_fits=export_fits,
                        )

