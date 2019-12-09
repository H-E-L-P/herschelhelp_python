import numpy as np
import matplotlib.pyplot as plt
import astropy
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from astropy.io import ascii, fits
from astropy.table import Table
from astropy import units as u

import .cutouts_server as ctts

def find_ra_and_dec(colnames):
    '''Finds the ra and dec column names regardless of what case they are in'''
    raname = ''
    decname = ''
    for n,colname in enumerate(colnames):
        if colname.lower() == 'ra':
            raname = colname
        if colname.lower() == 'dec':
            decname = colname
        if (raname != '') & (decname != ''):
            return(raname,decname)
        
def get_pix_length(image,wcs,box_length,ra_cent,dec_cent):
    #returns the number of pixels needed to encompass a box of given boxlength in degrees
    
    pix_centx, pix_centy = wcs.wcs_world2pix(ra_cent,dec_cent,0,ra_dec_order=True)
    pixscale_x, pixscale_y = np.array(wcs.wcs_pix2world(pix_centx,pix_centy,0,ra_dec_order=True)) - np.array(wcs.wcs_pix2world(pix_centx+1,pix_centy+1,0,ra_dec_order=True))

    pix_num_x = box_length/pixscale_x
    pix_num_y = box_length/pixscale_y
    
    return(abs(pix_num_x),abs(pix_num_y))

def cat_plot(ra,dec,wcs,ax,marker='dot',col='red',size=10,colz=False,redshift=None):
    '''
    Given a set of ra and dec will plot points on the given axis using the given wcs.
    The marker, size and colour of the points can be chosen and the points can be coloured
    according to their redshift
    -----------------------------
    Inputs
    ra: list or array of right ascensions
    dec: list or array of declinations
    wcs: the wcs of the axis the points are being displayed on
    ax: the matplotlib axis that the points are being displayed on
    marker: the marker style of the points (default dot/circle)
    col: The colour of the points (default red)
    size: The size of the points
    colz: whether to colour the points according to their redshifts
    redshift list or array of redshifts
    
    Outputs
    colbar: The colour bar corresponding to the redshifts so it can be plotted 
    
    '''
    pixx,pixy = wcs.wcs_world2pix(ra,dec,0,ra_dec_order=True)
    if colz==True:
        print('using redshift')
        col = redshift

    colbar = ax.scatter(ra,dec,marker=marker,c=col,s=size,cmap='Spectral',transform=ax.get_transform('world'))
    if colz==True:
        return(colbar)


def image_plot(image,ra,dec,wcs,box_length,fig):
    '''
    create a matplotlib image plot of a cutout of the given image
    --------------------------
    Inputs
    image: 2d array contianing the image data
    ra: right ascension of the center of the image
    dec: declination of the center of the image
    wcs: the wcs of the image
    box_length: give the boxlength in arcseconds
    fig: the figure that the matplotlib axis will be created in
    
    Output
    wcscut: the wcs of the cutout
    ax: the axis containg the cutout image
    '''
    cmap = 'binary'
    origin = 'lower'
    
    c = SkyCoord(ra*u.degree,dec*u.degree,unit='deg')
    length_x, length_y = get_pix_length(image,wcs,box_length,ra_cent,dec_cent)
    pix_cent_x, pix_cent_y = wcs.wcs_world2pix(ra,dec,0,ra_dec_order=True)
    pixscale = np.array(wcs.wcs_pix2world(0,0,0)) - np.array(wcs.wcs_pix2world(1,1,0))
    imgcut = Cutout2D(image, c, size=[box_length*u.degree,box_length*u.degree], wcs=wcs)
    wcscut = imgcut.wcs
    imgcut = imgcut.data

    vmin, vmax = ctts.cutout_scale(imgcut)
    ax = fig.add_subplot(111, projection=wcscut)

    ax.imshow(imgcut, vmax=vmax, vmin=vmin, cmap=cmap, interpolation='nearest', origin=origin)

    return(wcscut,ax)
    

def contour_plot(data,ra,dec,box_length,ax):
    '''
    create a contour plot of the given data on the given axis
    -------------------------
    Inputs
    data: a 2d array of the data that will be used to create the contour plot
    ra: the right ascension of the center of the data
    dec: the declination of the center of the data
    box_length: the length of the cutout in degrees
    ax: a matplotlib axis that the contour plot will be created in
    
    Output
    ax: the axis that now contains the contour plot
    '''
    
    image,wcs = data
    c = SkyCoord(ra*u.degree,dec*u.degree,unit='deg')
    length_x, length_y = get_pix_length(image,wcs,box_length,ra_cent,dec_cent)

    imgcut = Cutout2D(image, c, size=[box_length*u.degree,box_length*u.degree], wcs=wcs)
    wcscut = imgcut.wcs.celestial
    imgcut = imgcut.data
    vmin, vmax = ctts.cutout_scale(imgcut)
    print('contours are at {}'.format(np.linspace(vmin,vmax,7)))
    
    ax.contour(imgcut,levels=np.linspace(vmin,vmax,7),colors='white',transform=ax.get_transform(wcscut))
    
    return(ax)

def plot_figure(image,cat,wcs,ra_cent,dec_cent,plot_params,contour_data=None,return_fig=False):
    '''
    creates a figure containing a cutout of a image with points, from multiple catalogues, 
    overplotted on it which can be coloured according to their redshift. In addition a contour 
    plot of a sperate image can be overplotted
    --------------------
    Input
    image: A 2d array contaiing teh image data
    cat: an astropy table or list of astropy tables which contain the sources that you want 
    overplotted on the image. These catalogues do not need to be cutdown to the image area beforehand
    wcs: the wcs of the image being plotted
    ra_cent: the right ascension of the center of the desired cutout
    dec_cent: the declination of the center of the desired cutout
    plot_params: A dictionary containing the parameters for plotting the catalogues
        marker: a list of matplotlib markers for each catalogue
        col: a list of colors, one for each catalogue
        size: a list of sizes, one size for each catalogue
        box_length: the side legth of the desired cutout in degrees
        use_redshift: boolean for whether colours should be assigned based on the points redshift
    contour_data: default None. If not none should be a list containing [data,wcs] where data is a 2d array 
    containing the data that will be used to onstruct the contour plot and wcs is the wcs of the contour data
        
    '''
    
    #unpack plot_params
    marker,col,size,box_length,use_redshift = plot_params.values()
    
    #plot the cutout of the image and return the axis and the wcs of the cutout
    fig = plt.figure(figsize=(10,10))
    length_x, length_y = get_pix_length(image,wcs,box_length,ra_cent,dec_cent)
    wcsimg,ax = image_plot(image,ra_cent,dec_cent,wcs,box_length,fig)
    
    #create the radio contour plot and superimpose it onot the image cutout
    if contour_data != None:
        ax = contour_plot(contour_data,ra_cent,dec_cent,box_length,ax)
        ax.set_xlim(0,length_y-1)
        ax.set_ylim(0,length_y-1)
    
    #select the sources within the cutout from the catalogue and plot them
    ramax = ra_cent + box_length/2/np.cos(box_length/2)
    ramin = ra_cent - box_length/2/np.cos(box_length/2)

    decmax = dec_cent + box_length/2
    decmin = dec_cent - box_length/2
    #print('ramax = {}\nramin = {}\ndecmax = {}\ndecmin = {}'.format(ramax,ramin,decmax,decmin))
    if type(cat) != list:
        raname,decname = find_ra_and_dec(cat.colnames)
        mask = (cat[raname]>ramin) & (cat[raname]<ramax)  &  (cat[decname]>decmin) & (cat[decname]<decmax)
        cat = cat[mask]
        print('plotting {} sources'.format(len(cat)))
        if plot_params['use_redshift']==False:
            cat_plot(cat[raname],cat[decname],wcsimg,ax,marker[0],col,size)
        else:
            #the redshift column may need to be changed if your redhsift column isn't called 'redshift'
            colbar = cat_plot(cat[raname],cat[decname],wcsimg,ax,marker[0],col,size,colz=True,redshift=catalogue['redshift'])
            fig.colorbar(colbar)
    else:
        for n,catalogue in enumerate(cat):
            raname,decname = find_ra_and_dec(catalogue.colnames)
            mask = (catalogue[raname]>ramin) & (catalogue[raname]<ramax)  &  (catalogue[decname]>decmin) & (catalogue[decname]<decmax)
            catalogue = catalogue[mask]
            print('plotting {} sources'.format(len(catalogue)))
            if plot_params['use_redshift'][n]==False:
                cat_plot(catalogue[raname],catalogue[decname],wcsimg,ax,marker[n],col[n],size)
            else:
                #the redshift column may need to be changed if your redhsift column isn't called 'redshift'
                colbar = cat_plot(catalogue[raname],catalogue[decname],wcsimg,ax,marker[n],col[n],size,colz=True,redshift=catalogue['redshift'])
                fig.colorbar(colbar)
            
        
    if return_fig==True:
        return(fig)

    plt.show()
    