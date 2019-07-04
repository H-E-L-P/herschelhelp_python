# -*- coding: utf-8 -*-

import logging
import os

import numpy as np

import pyvo as vo

import yaml

from pymoc import MOC

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

from herschelhelp_internal.utils import inMoc
    
def clean_table(table):
    """Take a table produced by a VO query and remove all empty columns
    
    Often many columns are empty and make the tables hard to read.
    The function also converts columsn that are objects to strings.
    Object columns prevent writing to fits.
    
    Inputs
    =======
    table,    Astropy.table.Table
        The input table
    
    Returns
    =======
    table,    Astropy.table.Table
         The modified table.
    
    """
    table = table.copy()
    for col in table.colnames:
        #Remove empty columns
        if np.all(table[col].mask):
            print("Removing empty column: {}".format(col))
            table.remove_column(col)
            
        #Get rid of column type object from VO queries
        if table[col].dtype == 'object':
            print("Converting column {} type from object to string".format(col) )
            table[col] = table[col].astype(str)
 
        #Get rid of unit '-' from some tables
        if table[col].unit == '-':
            print("Converting column {} unit from '-' to None".format(col) )
            table[col].unit = None   
            
        #replace masked floats with nans     
        if table[col].dtype == float:
            table[col].fill_value = np.nan
    
    table = table.filled()
            
    return table
    
    
    
def help_cut_out(original, ra, dec, target=None, size_angle=100*u.arcsec):
    """Take an input image and cutout a square around the centre and 
    save the new file to the target
    
    
    Inputs
    ======
    original, str
        The location of the original image. It expects a HELP homogenised
        SPIRE, PACS map with the image in hdu 1, nebulised image in hdu 2
        and noise in hdu 3. For HS82 you may have to use the 0.9 map which 
        has no noise.
        
    ra, numpy.array
        An array of the right ascension
        
    dec, numpy.array
        An array of the declinations must be same length as ra
        
    target, str
        Where to write the output. If None just returns the image without writing.
    
    
    
    """
    original_fits = fits.open(original)
    primary_hdu = original_fits[0]
    image_hdu = original_fits[1]
    noise_hdu = original_fits[2]
    image_wcs = WCS(image_hdu.header)
    noise_wcs = WCS(noise_hdu.header)
    
   # Make the cutout, including the WCS
    ra_dec_position = SkyCoord(ra, dec, frame='icrs')
    x_y_position = skycoord_to_pixel(ra_dec_position, image_wcs)
    
    #Some of images use CD2_2 key, some use CDELT2
    try:
        pix_size_deg = original_fits[1].header['CD2_2']
    except KeyError:
        pix_size_deg = original_fits[1].header['CDELT2']
        
    size_pix = (
        2*(Angle(size_angle).degree/
        Angle(pix_size_deg*u.deg ).degree )*  u.pixel
        )
    image_cutout = Cutout2D(
        image_hdu.data, 
        position=x_y_position, 
        size=size_pix, 
        wcs=image_wcs
        )
    noise_cutout = Cutout2D(
        noise_hdu.data, 
        position=x_y_position, 
        size=size_pix, 
        wcs=noise_wcs
        )

    # Put the cutout image in the FITS HDU
    image_hdu.data = image_cutout.data
    noise_hdu.data = noise_cutout.data

    # Update the FITS header with the cutout WCS
    image_hdu.header.update(image_cutout.wcs.to_header())
    noise_hdu.header.update(noise_cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    
    hdu_list = fits.HDUList([primary_hdu, image_hdu, noise_hdu])
    
    if target != None:
        hdu_list.writeto(target, overwrite=True)
    hdu_list.close()
    
    return target, image_hdu.data, noise_hdu.data
    
def help_id_to_ra_dec(help_id):
    """Take a HELP id and return the ra dec in deg
    
    TODO: vectorise
    
    parameters
    ----------
    
    help_id: string
       The help ID
       
    returns
    -------
    
    ra: float
       Right Ascension in degrees
       
    dec: float
        Declination in degrees
    """
    if '-' in help_id:
        ra, dec = help_id.strip('HELP_J').split('-')
        dec = '-' + dec
    if '+' in help_id:
        ra, dec = help_id.strip('HELP_J').split('+')
        dec = '+' + dec
    
    ra = ra[:2] + 'h' + ra[2:]
    ra = ra[:5] + 'm' + ra[5:]
    ra = ra + 's' 
    
    
    
    dec = dec[:3] + 'd' + dec[3:]
    dec = dec[:6] + 'm' + dec[6:]
    dec = dec + 's'
    
    coords = SkyCoord('{} {}'.format(ra, dec))
    
    return coords.ra.value, coords.dec.value
    
    
def position_to_field(ra, dec, dmu_loc='../dmu_products/'):
    """Take a position and return HELP field names
    
    If  not in HELP return None. This does not check if the
    ID is in DR1 or a later DR. It simply converts the ID to ra
    dec values and checks if they are in a field moc
    
    
    Inputs
    ======
    
    ra, np.array
        Array of right ascension positions
        
    dec, np.array
        Array of declination positions
        
    Returns
    =======
    
    field, str
        The HELP field. 'not_help' if outside HELP.
    """
    field = None
    try:
        fields_info = yaml.load(open(dmu_loc + "dmu2/meta_main.yml", 'r'))
    except FileNotFoundError:
        print("{} is not present or does not contain the"
              + " field description YML file.")
        return None
 
    field = np.full(len(ra), 'not_help', dtype='<U18')
    for f in fields_info['fields']:
        f_moc = MOC(filename=dmu_loc + 
                    "dmu2/dmu2_field_coverages/{}_MOC.fits".format(f['name']))
        in_field = inMoc(ra, dec, f_moc)
     
        field[in_field] = f['name']

    return field