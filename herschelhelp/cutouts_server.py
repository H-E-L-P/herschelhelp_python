"""
    Author: Estelle Pons epons@ast.cam.ac.uk
    Date: 2018
"""

# Python 2/3 compatibility
from __future__ import print_function   # to use print() as a function in Python 2

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import os
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u
import scipy.ndimage as ndimage
from astropy.wcs import WCS
import astropy.io.fits as fits
from astropy.table import Table

import cutouts_dl as cdl


##############################################################################
##                            Read & Plot Functions                         ##
##############################################################################
def radecStr(ra, dec, precision=1):
    """
    make a radec name from ra, dec e.g. HHMMSSsDDMMSS  
    """   
    radec = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    radecname = radec.to_string('hmsdms', decimal=False, sep='', precision=precision)
    radec_str = np.core.defchararray.replace(radecname, ' ', '')
    radec_str = str(radec_str)
        
    return radec_str
        
        
def cutout_scale(im, num_min = 2.0, num_max = 5.0):
    """
    Takes an image array and returns the vmin and vmax required to scale the image 
    between median + 5 * sigma MAD and median - 2 * sigma MAD
    """

    import numpy as np
    import astropy.stats as apys

    data = im.flatten()   #Return a copy of the array collapsed into one dimension

    try:
        med = np.median(data[np.isnan(data)!=1])
        sigma_MAD = apys.mad_std(data[np.isnan(data)!=1])
    except IndexError:
        med = 0.0
        sigma_MAD = 0.0
    vmax = med + num_max * sigma_MAD
    vmin = med - num_min * sigma_MAD

    return vmin, vmax
    

def urlfile_exists(location):
    try:
        import urllib2
        httpErr = urllib2.HTTPError
    except:
        import urllib
        import urllib.request as urllib2
        httpErr = urllib.error.HTTPError

    if location != "":
        try:
            ret = urllib2.urlopen(location)
            return True
        except httpErr:
             return False
        except ValueError:
             return False

   
def rd_fits(filename, ra, dec, hdrNum=1, width_as=20., pixelscale=None, hdrKey_pixelscale=None, smooth=False):
    """
    Read fits image and create cutout of size width_as
    """
    ## read fits
    if os.path.exists(filename) or urlfile_exists(filename):
        print("     fits file available")
        h = fits.open(filename, cache=False)
        if hdrNum is not None:
            hdr = h[hdrNum].header 
            data = h[hdrNum].data
        else:
            # CFHT multi-extension files
            print("     len(hdulist):", len(h))
            for hnum in range(len(h)):
                print("     hdrNum:", hnum)
                data_hnum = h[hnum].data
                if data_hnum is not None and np.nanmean(data_hnum) !=0:
                    hdr = h[hnum].header
                    data = data_hnum
                    if "NAXIS3" in hdr:
                        data = data[0]          
                    break
        h.close()
            
        try:
            w = WCS(hdr, naxis = 2)
        except:    
            # HST fits: crop header (to prevent issue with WCS)
            if "D001VER" in hdr.keys():
                new_keys = ["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "TELESCOP", "INSTRUME", "EQUINOX",\
                            "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"]
                new_hdul = fits.PrimaryHDU()
                for new_key in new_keys:
                    new_hdul.header[new_key] = hdr[new_key]
                hdr = new_hdul.header
                w = WCS(hdr, naxis = 2)
                             
        pix_coord = w.wcs_world2pix(ra, dec, 1)


        if pixelscale is None:
            pixelscale = hdr[hdrKey_pixelscale]
        #print(pixelscale
        width_pix = width_as / pixelscale

        cutout = Cutout2D(data, pix_coord, (width_pix, width_pix))
        image = cutout.data

        #Smooth the data
        if smooth:
            image = ndimage.gaussian_filter(image, sigma=1.0, order=0)
            
    else:
        print("     No fits file available")
        null_image = np.zeros(shape=(10,10), dtype="int8")
        image = null_image

    return image
    
    
def plt_image(band, image, fig, ax, psfmags=None, cmap="binary", minmax="MAD", origin="lower"):
    """
    Plot the image (cutout size defined in rd_fits)
    """
    ## Plot min - max
    if minmax == "MAD":
        vmin, vmax = cutout_scale(image)
    else:
        vmin, vmax = np.nanmin(image), np.nanmax(image)

    ## imshow
    ax.imshow(image, vmax=vmax, vmin=vmin, cmap=cmap, interpolation='nearest', origin=origin)
    
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    
    if psfmags is not None:
        ax.set_title('{} (mag = {:.2})'.format(band, psfmags))
    else:
        ax.set_title(band)
    
    return vmin, vmax
    
    
##############################################################################
##                              Surveys Functions                           ##
##############################################################################

# -------------------------------------------------------------------------- #
#                                  2MASS                                     #
# -------------------------------------------------------------------------- #
def cutout_twomass(ra, dec, bands=['J','H','K'], psfmags=None,\
                  imDir="/data/2mass/", input_filename=[], getFullIm=False, saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                  return_val=False, saveDir=None):

    """
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image color map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: "2MASS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## 2MASS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: filename of the fits file if save on disk
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("2MASS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.twomass_dl(ra, dec, band, width_as=width_as, FitsOutputPath=imDir,\
                                    getFullIm=getFullIm, saveFITS=saveFITS)
                                    
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=1., smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
            
        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))
        
    ## Add a title to the figure
    if figTitle:
        fig.suptitle('2MASS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                     ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None: 
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_2MASS_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig


# -------------------------------------------------------------------------- #
#                                   CFHT                                     #
# -------------------------------------------------------------------------- #
def cutout_cfht(ra, dec, bands=["u", "g", "r", "i", "z"], instrument="MegaPrime", optFiltersGen='both', psfmags=None,\
                imDir="/data/cfht/", input_filename=[], getFullIm=False, saveFITS=False,\
                width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                return_val=False, saveDir=None):

    """
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image color map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: CFHT cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## CFHT parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts (ugriz for MegaPrime; YJHK for WIRCam; BRI for CFHT12k)
        instrument: CFHT intrument (MegaPrime, WIRCam, CFH12K MOSAIC)
        optFiltersGen: version of the MegaPrime filters (1stGen, 3rdGen or both)
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: filename of the fits file if save on disk
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("CFHT cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.cfht_dl(ra, dec, band, instrument=instrument, width_as=width_as, FitsOutputPath=imDir,\
                                    getFullIm=getFullIm, saveFITS=saveFITS)
                                    
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=None, width_as=width_as, hdrKey_pixelscale="PIXSCAL1", smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
            
        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))
        
    ## Add a title to the figure
    if figTitle:
        fig.suptitle('CFHT cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                     ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None: 
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_CFHT_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
        
        
# -------------------------------------------------------------------------- #
#                                  DECaLS                                    #
# -------------------------------------------------------------------------- #
def cutout_decals(ra, dec, bands=['g','r','z'], dr=5, psfmags=None, getFullIm=False, brickid=None,\
                  imDir="/data/vault/des/DECALS/DR5/coadd/", saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                  return_val=False, saveDir=None):

    """
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image color map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: DECaLS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## DECaLS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        dr: DECaLS data release
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        brickid: brickid of the DECaLS observation, i.e. in the form "2054p095"
        getFullIm: get the full DECaLS tile if True
        imDir: directory of the fits file if already save on disk
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("DECaLS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        brick = str(brickid)[0:3]
        filename = imDir + str(brick) + "/" + str(brickid) + "/" + "legacysurvey-{}-image-{}.fits.fz".format(str(brickid), band)
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename):
            filename, brickid = cdl.decals_dl(ra, dec, band, dr=dr, width_as=width_as, getFullIm=getFullIm,\
                                       brickid=brickid, FitsOutputPath=imDir, saveFITS=saveFITS)
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        if getFullIm:
            hdrNum = 1
        else:
            hdrNum = 0
        
        image = rd_fits(filename, ra, dec, hdrNum=hdrNum, width_as=width_as, pixelscale=0.262, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
            
        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))
        
    ## Add a title to the figure
    if figTitle:
        fig.suptitle('DECalS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                     ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None: 
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_DECaLS_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
        

# -------------------------------------------------------------------------- #
#                                    DES                                     #
# -------------------------------------------------------------------------- #
def cutout_des(ra, dec, bands=['g','r','i','z','y'], dr=1, psfmags=None,\
                  imDir="/des/", input_filename=[], saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                  return_val=False, saveDir=None):
    """
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image color map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: DES cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## DES parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        dr: DES data release
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: names of the fits filename for each band (if save on disk) 
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """ 
    
    print("DES-DR{:d} cutout(s), band(s):".format(dr), "".join(bands))              
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
   
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        pixelscale = 0.27
        if not os.path.exists(filename):
            filename, pixelscale = cdl.des_dl(ra, dec, band, dr=dr, width_as=width_as, FitsOutputPath=imDir, saveFITS=saveFITS)

        print("   filename:", filename)
        
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=pixelscale, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('DES-DR{:d} cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(dr, width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_DES-DR{:d}_{}_{}_{:.0f}arcsec.png".format(dr, radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
                      

# -------------------------------------------------------------------------- #
#                                 FLS - KPNO                                 #
# -------------------------------------------------------------------------- # 
def cutout_flsKPNO(ra, dec, bands=['R'], psfmags=None, imDir="/data/fls-kpno/", input_filename=[], saveFITS=False,\
               width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True,\
               return_val=False, saveDir=None):   
    """
        Plot the FLS-KPNO R-band cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: FLS-KPNO cutout 20"x20" ra=, dec= (Jradec); default is True)
        
        ## FLS-KPNO parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: FLS-KPNO bands (only R)
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: names of the fits filename for each band (if save on disk) 
        saveFITS: save fits file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("FLS-KPNO cutout, R-band:")
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
   
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename):
            filename = cdl.flsKPNO_dl(ra, dec, width_as=width_as, FitsOutputPath=imDir, saveFITS=saveFITS)

        print("   filename:", filename)
        
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, hdrKey_pixelscale="PIXSCAL1", smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('FLS-KPNO cutout ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_FLS-KPNO_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
        
        
        
# -------------------------------------------------------------------------- #
#                     Hawaii Hubble Deep Field North                         #
# -------------------------------------------------------------------------- #
def cutout_hdfn(ra, dec, bands=["U","B","V","R","I","z","HK"], psfmags=None, 
                  imDir="/data/hdfn/", saveFITS=False, getFullIm=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True,\
                  return_val=False, saveDir=None):
    """
        Plot all the HDFN bands cutouts on one plot for an input source position
        http://www.astro.caltech.edu/~capak/hdf/
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: PanSTARRS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## HDFN parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        saveFITS: save fits file on disk (to imDir)
        getFullIm: get the full fits image if True
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """

    print("HDFN cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
   
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if band == "V":
            band == "V0201"
            
        if getFullIm:
            filename = imDir + "HDF.{}.fits.gz".format(band)
        else:
            filename = imDir + "{}-{}_{}_{:.0f}arcsec.fits".format(radec_str, "HDFN", band, width_as)
        
        if ra > 188.74 and ra < 189.70 and dec > 61.96 and dec < 62.41:
            ### If filename does nor exists -> get file from url
            if not os.path.exists(filename):
                filename = cdl.hdfn_dl(ra, dec, band, width_as=width_as, FitsOutputPath="/data/hdfn/",\
                                       getFullIm=getFullIm, saveFITS=saveFITS)

            print("   ", filename)
        else:
            filename = ""
            print("NO HDFN {}-BAND COVERAGE FOR THIS POSITION".format(band))
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=0.3, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('HDFN cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_HDFN_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
        


# -------------------------------------------------------------------------- #
#                                   HSC                                      #
# -------------------------------------------------------------------------- #
def cutout_hsc(ra, dec, bands=["g","r","i","z","y","N816","N921"], login=[],\
                psfmags=None, imDir="/data/hsc/", input_filename=[], saveFITS=False,\
                width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True,\
                return_val=False, saveDir=None):
    """
        Plot all the HSC bands cutouts on one plot for an input source position
        An account is required to access public data
          SEE https://hsc-release.mtk.nao.ac.jp/doc/index.php/tools/
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: HSC cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## HSC parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        login: username and password array. COMPULSORY to access HSC public dataset
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: names of the fits filename for each band (if save on disk) 
        saveFITS: save fits file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """

    print("HSC cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
   
    if len(login) == 0:
        print("ACCOUNT REQUIRED TO ACCESS HSC PUBLIC DATA (see https://hsc-release.mtk.nao.ac.jp/doc/index.php/tools/)")
        return
        
    else:
        for i, band in enumerate(bands):
            print("{}-band".format(band))
            ### Filename of fits image if save of the disk
            if len(input_filename) == 0: 
                input_filename = ""
            else:
                input_filename = input_filename[i]
            filename = imDir + input_filename
        
            ### If filename does nor exists -> get file from url
            if not os.path.exists(filename):
                filename = cdl.hsc_dl(ra, dec, band, width_as=width_as, login=login, FitsOutputPath=imDir, saveFITS=saveFITS)

            print("   ", filename)

            ### Read fits file: cutout size = width_as
            ###                 filename could be a system path or an url or ""
            print("   Try to read the fits file ...")
            
            image = rd_fits(filename, ra, dec, hdrNum=1, width_as=width_as, pixelscale=0.17, smooth=smooth)
            
            ### Plot image: cutout size = width_as
            print("   Plot the cutout ...")
            ax = fig.add_subplot(gs[0,i])
            
            if psfmags is not None:
                psfmags = psfmags[i]

            vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
            datas.append((image, vmin, vmax))

        ## Add a title to the figure
        if figTitle:
            fig.suptitle('HDFN cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                                                ra, dec, radec_str), fontsize=15)
            
        ### Output
        if return_val:
            print(" Return image data")
            plt.close(fig)
            return datas
        
        if saveDir is not None:
            print(" Save the figure to", saveDir)
            allBands = "".join(bands)
            plt.savefig(saveDir + "Cutouts_HDFN_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                        bbox_inches="tight")
            plt.close()
        else:
            print(" Return the figure")
            return fig
        
        
# -------------------------------------------------------------------------- #
#                             HST/ACS - HST/WFC3                             #
# -------------------------------------------------------------------------- # 
def cutout_hst(ra, dec, bands=["F435W","F606W","F775W","F814W","F850LP"], instrument="ACS",\
               psfmags=None, imDir="/data/hst/", input_filename=[], saveFITS=False,\
               width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True,\
               return_val=False, saveDir=None):   
    """
        Plot all the HST/ACS or HST/WFC3 bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: PanSTARRS cutout 20"x20" ra=, dec= (Jradec); default is True)
        
        ## HST parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        instrument: HST instrument ("ACS" = optical or "WFC3" = IR)
            ["F435W","F606W","F775W","F814W","F850LP"] for ACS (B Johnson, V WFPC2, i SDSS, I WFPC2, z SDSS respectively)
            ["F105W","F125W","F140W","160W","F098M"] for WFC3 (Wide Z, J, JH, Hs, blue grusm ref respectively)
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: names of the fits filename for each band (if save on disk) 
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        saveFITS: save fits file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("HST/{} cutout(s), band(s):".format(instrument), "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
   
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename):
            filename = cdl.hst_dl(ra, dec, band, instrument=instrument, width_as=width_as,\
                                  FitsOutputPath=imDir, saveFITS=saveFITS)

        print("   filename:", filename)
        
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        pixelscales = {"ACS": 0.05, "WFC3": 0.09}
        pixelscale = pixelscales[instrument]

        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=pixelscale, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('HST/{} cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(instrument, width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_HST-{}_{}_{}_{:.0f}arcsec.png".format(instrument, radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
    
    
# -------------------------------------------------------------------------- #
#                              Legacy Survey                                 #
# -------------------------------------------------------------------------- #
def cutout_legacySurvey(ra, dec, bands=['g','r','z'], dr=6, psfmags=None,\
                  imDir="/data/LegacySurvey/", input_filename=[], saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                  return_val=False, saveDir=None):

    """
        Legacy Survey = BASS (gr-bands) + MzLS (z-band)
        
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image color map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: Legacy Survey cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## Legacy Survey parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        dr: Legacy Survey data release
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("Legacy Survey cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename):
            filename = cdl.legacySurvey_dl(ra, dec, band, dr=dr, width_as=width_as,\
                                           FitsOutputPath=imDir, saveFITS=saveFITS)
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=0.262, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
            
        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))
        
    ## Add a title to the figure
    if figTitle:
        fig.suptitle('Legacy Survey (BASS+MzLS) cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                     ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None: 
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_LegacySurvey_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
        


# -------------------------------------------------------------------------- #
#                                  NDWFS                                     #
# -------------------------------------------------------------------------- #
def cutout_ndwfs(ra, dec, bands=['Bw','R','I','K'], psfmags=None,\
                  imDir="/data/ndwfs/", input_filename=[], saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                  return_val=False, saveDir=None):

    """
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image color map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: "NDWFS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## NDWFS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: filename of the fits file if save on disk
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("NDWFS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.ndwfs_dl(ra, dec, band, width_as=width_as, FitsOutputPath=imDir, saveFITS=saveFITS)
                                    
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=None, hdrKey_pixelscale="PIXSCAL1", smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
            
        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))
        
    ## Add a title to the figure
    if figTitle:
        fig.suptitle('NDWFS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                     ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None: 
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_NDWFS_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig


            
# -------------------------------------------------------------------------- #
#                               Pan-STARRS                                   #
# -------------------------------------------------------------------------- #
def cutout_ps1(ra, dec, bands=["g","r","i","z","y"], psfmags=None, 
                  imDir="/data/ps1ardata/ps1dr1_cutouts_cache/", getFullIm=False, saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True,\
                  return_val=False, saveDir=None):
    """
        Plot all the PanSTARRS bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: PanSTARRS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## Pan-STARRS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        saveFITS: save fits file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """

    print("PanSTARRS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.95, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    ### Filename of fits image if save of the disk (20" cutouts)
    dic_list = {band+"im": glob.glob(full_imDir + radec_str + '*.unconv.fits'.format(band)) for band in bands}

    projcell, subcell = None, None
    
    bandIms = [band+"im" for band in bands]
    for (i, bandIm) in enumerate(bandIms):
        print("{}-band".format(bands[i]))
        
        ### If filename does nor exists -> get file from url
        if len(dic_list[bandIm]) == 0 or width_as > 20.:
            filename, projcell, subcell = cdl.ps1_dl(ra, dec, bands[i], radec_str, projcell=projcell, subcell=subcell,\
                                                     getFullIm=getFullIm, width_as=width_as, saveFITS=saveFITS,\
                                                     FitsOutputPath=imDir)
        else:
            filename = dic_list[bandIm][0]
         
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=0.258, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
            
        vmin, vmax = plt_image(bands[i], image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))
        
    ## Add a title to the figure
    if figTitle:
        fig.suptitle('PanSTARRS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                     ra, dec, radec_str), fontsize=15)
                     
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None: 
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_PS1_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig

# -------------------------------------------------------------------------- #
#                                  SDSS                                      #
# -------------------------------------------------------------------------- # 
def cutout_sdss(ra, dec, bands=['u','g','r','i','z'], dr=12, objid=None, psfmags=None,\
                imDir="/data/sdss/", input_filename=[], saveFITS=False,\
                width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                return_val=False, saveDir=None):
    """
        Plot all the SDSS bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: SDSS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## SDSS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        dr: SDSS data release (default is 12)
        objid: SDSS objid 
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        saveFITS: save fits file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """                   
            
    print("SDSS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
    params = {}
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            if len(params) == 0:
                filename, params = cdl.sdss_dl(ra, dec, band, dr=dr, objid=objid, FitsOutputPath=imDir, saveFITS=saveFITS)
            else:
                filename = "{base_url}frame-{}-{run:06d}-{camcol}-{field:04d}.fits.bz2".format(band, **params)

        print("   ", filename)
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, pixelscale=0.396, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('SDSS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_SDSS_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig
        
            
            
            
                   
# -------------------------------------------------------------------------- #
#                                  Spitzer                                   #
# -------------------------------------------------------------------------- # 
def cutout_spitzer(ra, dec, bands=['I1','I2','I3','I4'], dataset="SEIP", psfmags=None,\
                   imDir="/data/spitzer/", input_filename=[], saveFITS=False,\
                   width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                   return_val=False, saveDir=None):
    """
        Plot all the Spitzer/IRAC bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: Spitzer cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## Spitzer parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        dataset: Spitzer dataset -> SEIP (default; includes SDWFS, SWIRE, SERVS), SDWFS, SERVS, SHELA, SIMES, SSDF, SpIES, SpUDS, SWIRE 
        imDir: directory of the fits file if already save on disk
        saveFITS: save fits file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
                    
    print("Spitzer-{}/IRAC cutout(s), band(s):".format(dataset), "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.spitzer_dl(ra, dec, band, width_as=width_as, dataset=dataset,\
                                                FitsOutputPath=imDir, saveFITS=saveFITS)

        print("   ", filename)
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        ### pixelscale mosaic (SEIP, SSDF, SpIES, SpUDS) = 0.6
        print("   Try to read the fits file ...")
        pixelscales = {"SEIP":None, "SDWFS":0.84, "SERVS":0.6, "SHELA":0.8, "SSDF":None, "SpIES":None, "SpUDS":None,\
                       "SWIRE":0.6, "GOODS":0.6, "FLS":1.22}
        hdrKey_pixelscales = {"SEIP": "PIXSCALE", "SDWFS": None, "SERVS": None, "SHELA": None, "SSDF": "PXSCAL2",\
                              "SpIES": "PXSCAL2", "SpUDS": "PXSCAL2", "SWIRE": None, "GOODS": None, "FLS": None}
        image = rd_fits(filename, ra, dec, hdrNum=0, width_as=width_as, smooth=smooth,\
                        pixelscale=pixelscales[dataset], hdrKey_pixelscale=hdrKey_pixelscales[dataset])
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('Spitzer-{}/IRAC cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(dataset, width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_Spitzer-{}_{}_{}_{:.0f}arcsec.png".format(dataset, radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig          
                   
 
# -------------------------------------------------------------------------- #
#                                    UHS                                     #
# -------------------------------------------------------------------------- #
def cutout_uhs(ra, dec, bands=['J'], database="UHSDR1", wsaLogin=[],\
               psfmags=None, imDir="/data/uhs/", input_filename=[], saveFITS=False,\
               width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
               return_val=False, saveDir=None):
    
    """
        Plot all the bands cutouts on one plot for an input source position
        
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: UHS cutout 20"x20" ra=, dec= (Jradec); default is True)
        
        ## UHS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                 Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: name of the input file if save on disk
        database: UHS database used "UHSDR{}" (default if DR1)
        wsaLogin: array of [username, password] (Not compulsory, data puplic since 01/08/2018)
        saveFITS: save fits tile file on disk (to imDir)
        
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("UHS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.uhs_dl(ra, dec, band, wsaLogin=wsaLogin, database=database, width_as=width_as,\
                                  FitsOutputPath=imDir, saveFITS=saveFITS)

        print("   ", filename)
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=1, width_as=width_as, hdrKey_pixelscale="PIXLSIZE", smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('UHS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_UHS_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig


# -------------------------------------------------------------------------- #
#                                  UKIDSS                                    #
# -------------------------------------------------------------------------- #
def cutout_ukidss(ra, dec, bands=['Z','Y','J','H','K'], programme="LAS", database="UKIDSSDR10PLUS", wsaLogin=[], 
                  psfmags=None, imDir="/data/ukidss/", input_filename=[], saveFITS=False,\
                  width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                  return_val=False, saveDir=None):

    """
        Plot all the bands cutouts on one plot for an input source position
    
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: UKIDSS cutout 20"x20" ra=, dec= (Jradec); default is True)
    
        ## UKIDSS parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: name of the input file if save on disk
        programme: survey name, i.e DXS, GCS, GPS, LAS, UDS (default is LAS)
        database: UKIDSS database used "UKIDSSDR{}PLUS" (default if DR10)
        wsaLogin: array of [username, password, community] to access non-public data
        saveFITS: save fits tile file on disk (to imDir)
    
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("UKIDSS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))

    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.ukidss_dl(ra, dec, band, wsaLogin=wsaLogin, programme=programme, database=database,\
                                     width_as=width_as, FitsOutputPath=imDir, saveFITS=saveFITS)
                
        print("   ", filename)
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        ### pixelscale = 0.4
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=1, width_as=width_as, hdrKey_pixelscale="PIXLSIZE", smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('UKIDSS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                        ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
                                                                                                    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_UKIDSS_{}_{}_{:.0f}arcsec.png".format(radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig



# -------------------------------------------------------------------------- #
#                                   VISTA                                    #
# -------------------------------------------------------------------------- #
def cutout_vista(ra, dec, bands=["Z","Y","J","H","K"], survey="VHS", database="VHSDR5", wsaLogin=[],\
               psfmags=None, imDir="/data/vista/", input_filename=[], saveFITS=False,\
               width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
               return_val=False, saveDir=None):
    
    """
        VISTA surveys: VHS, VIDEO, VIKING
        Plot all the bands cutouts on one plot for an input source position
        
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: VISTA cutout 20"x20" ra=, dec= (Jradec); default is True)
        
        ## VISTA parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts (VHS: JHK, VIKING: ZYJHK, VIDEO: ZYJHK
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
                 Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: name of the input file if save on disk
        survey: VISTA survey (VHS, VIDEO, VIKING)
        database: VISTA database used = survey + DataRealease (survey + data for proprietary data)
            Last public Release (up to June 2018): VHSDR5, VIKINGDR4, VIDEODR5
            full list: http://horus.roe.ac.uk:8080/vdfs/VgetImage_form.jsp
        wsaLogin: array of [username, password, community] to access non-public data
        saveFITS: save fits tile file on disk (to imDir)
        
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
    """
    
    print("VISTA-{} cutout(s), band(s):".format(survey), "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)

    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0: 
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.vista_dl(ra, dec, band, wsaLogin=wsaLogin, database=database, width_as=width_as,\
                                  FitsOutputPath=imDir, saveFITS=saveFITS)

        print("   ", filename)
    
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=1, width_as=width_as, pixelscale=0.34, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]

        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('VISTA-{} cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(survey, width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
        
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_VISTA-{}_{}_{}_{:.0f}arcsec.png".format(survey, radec_str, allBands, width_as),\
                    bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig


# -------------------------------------------------------------------------- #
#                                 VST-ATLAS                                  #
# -------------------------------------------------------------------------- #
def cutout_vstAtlas(ra, dec, bands=["u","g","r","i","z"], database="ATLASDR3",\
                 psfmags=None, imDir="/data/vst-atlas/", input_filename=[], saveFITS=False,\
                 width_as=20., smooth=False, cmap="binary", minmax="MAD", origin="lower", figTitle=True, \
                 return_val=False, saveDir=None):
    
    """
        Plot all the bands cutouts on one plot for an input source position
        
        ## Cutouts parameters
        width_as: size of the cutout box; default is 20arcsec
        smooth: gaussian smoothing with sigma=1.0; defaul is False
        cmap: image colour map
        minmax: Defined the min-max scale of the image; default is from sigma_MAD(image) (SEE def cutout_scale)
        origin: where to place the [0,0] index of the image; default is "lower"
        figTitle: add a title to the final figure (ex: VISTA cutout 20"x20" ra=, dec= (Jradec); default is True)
        
        ## VISTA parameters
        ra, dec: position of the source in deg (single object, not an array)
        bands: filters for which to do the cutouts
        psfmags: magnitudes of the source. Should be an array of the same size than bands or None (default)
        Will be added to band cutout title if not None
        imDir: directory of the fits file if already save on disk
        input_filename: name of the input file if save on disk
        database: ATLAS database used = ATLAS + DataRealease
        saveFITS: save fits tile file on disk (to imDir)
        
        ## Output parameters
        return_val: return image data, min-max(image); default is False
        saveDir: output directory to save the final figure. If None do not save; default is None
        """
    
    print("VST-ATLAS cutout(s), band(s):", "".join(bands))
    
    ### radec: HHMMSSsDDMMSS
    radec_str = radecStr(ra, dec, precision=1)
    
    ### Figure: defined fig and gs
    figWidth = len(bands) * 8./3.
    fig = plt.figure(figsize=(figWidth, 4))
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.90, bottom = 0, wspace = 0)
    gs = gridspec.GridSpec(1, len(bands))
    
    datas = []
    
    for i, band in enumerate(bands):
        print("{}-band".format(band))
        ### Filename of fits image if save of the disk
        if len(input_filename) == 0:
            input_filename = ""
        else:
            input_filename = input_filename[i]
        filename = imDir + input_filename
        
        ### If filename does nor exists -> get file from url
        if not os.path.exists(filename) or input_filename == "":
            filename = cdl.vstAtlas_dl(ra, dec, band, database=database, width_as=width_as,\
                                    FitsOutputPath=imDir, saveFITS=saveFITS)
        
        print("   ", filename)
        
        ### Read fits file: cutout size = width_as
        ###                 filename could be a system path or an url or ""
        print("   Try to read the fits file ...")
        image = rd_fits(filename, ra, dec, hdrNum=1, width_as=width_as, pixelscale=0.21, smooth=smooth)
        
        ### Plot image: cutout size = width_as
        print("   Plot the cutout ...")
        ax = fig.add_subplot(gs[0,i])
        
        if psfmags is not None:
            psfmags = psfmags[i]
        
        vmin, vmax = plt_image(band, image, fig, ax, psfmags=psfmags, cmap=cmap, minmax=minmax, origin=origin)
        datas.append((image, vmin, vmax))

    ## Add a title to the figure
    if figTitle:
        fig.suptitle('VST-ATLAS cutouts ({:.0f}"x{:.0f}") \n ra: {:.4f}, dec: {:.4f} (J{})'.format(width_as, width_as,\
                                                                                            ra, dec, radec_str), fontsize=15)
    
    ### Output
    if return_val:
        print(" Return image data")
        plt.close(fig)
        return datas
                                                                                                    
    if saveDir is not None:
        print(" Save the figure to", saveDir)
        allBands = "".join(bands)
        plt.savefig(saveDir + "Cutouts_VISTA-{}_{}_{}_{:.0f}arcsec.png".format(survey, radec_str, allBands, width_as),\
                        bbox_inches="tight")
        plt.close()
    else:
        print(" Return the figure")
        return fig

