# Python 2/3 compatibility
from __future__ import print_function   # to use print() as a function in Python 2

import os, sys
import urllib, re
try:
    import urllib2, urllib
    import cookielib
except:
    import urllib
    import urllib.request as urllib2
    import http.cookiejar as cookielib
import numpy as np
from astropy import table
from astropy.table import Table
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import xml.etree.ElementTree as ET



def radecStr(ra, dec, precision=1):
    """
        make a radec name from ra, dec e.g. HHMMSSsDDMMSS
        """
    radec = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    radecname = radec.to_string('hmsdms', decimal=False, sep='', precision=precision)
    radec_str = np.core.defchararray.replace(radecname, ' ', '')
    radec_str = str(radec_str)
    
    return radec_str


# -------------------------------------------------------------------------- #
#                                  2MASS                                     #
# -------------------------------------------------------------------------- #
def twomass_dl(ra, dec, band, getFullIm=False, width_as=20.,\
             FitsOutputPath="/data/2mass/", saveFITS=False):
    """
        Download 2MASS fits image
        https://irsa.ipac.caltech.edu/ibe/docs/twomass/allsky/allsky/
        
        Full image Located to :
            https://irsa.ipac.caltech.edu/ibe/data/twomass/allsky/allsky/{ordate:6d}{hemisphere:1s}/s{scanno:03d}/image/{fname:s} 
            ex: https://irsa.ipac.caltech.edu/ibe/data/twomass/allsky/allsky/990916n/s028/image/hi0280150.fits.gz
            
         ordate, hemisphere, scanno, and fname from:
             https://irsa.ipac.caltech.edu/ibe/search/twomass/allsky/allsky?POS{ra},{dec}
        
        Cutouts -> append center and size parameters at the end of the url
                   size units can be pixels (px, pix, pixels) or angular (arcsec, arcmin, deg, rad)
            ex:  https://irsa.ipac.caltech.edu//ibe/data/twomass/allsky/allsky/990916n/s028/image/hi0280150.fits.gz?center=352.38,38.63&size=100pix
            
       ------------------------------------------------------------------------------------------------------------------------
       
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        width_as: size of the cutout if getFullIm=False
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk (to FitsOutputPath/BRICKID[0:3]/BRICKID/) 
    """
    
    print("   Get the url of the fits file ...")
    ## Get  ordate, hemisphere, scanno, and fname
    url_search = "https://irsa.ipac.caltech.edu/ibe/search/twomass/allsky/allsky?POS=" + str(ra) + "," + str(dec)
    t = ascii.read(url_search)
    if len(t) > 0:
        t = t[0]
        params = {'ordate': t["ordate"], 'hemisphere': t["hemisphere"], 'scanno': t["scanno"],\
                  'fname': band.lower() + t["fname"][1:]}
    else:
        print("NO 2MASS {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return ""
        
    ## Get the fits image
    url = "https://irsa.ipac.caltech.edu/ibe/data/twomass/allsky/allsky/"
    path = '{ordate:6s}{hemisphere:1s}/s{scanno:03d}/image/{fname:s}'.format(**params)
    urlIm = url + path

    if not getFullIm:
        # Get a cutout of size width_as
        urlIm = urlIm + "?center={:f},{:f}&size={:.0f}arcsec".format(ra, dec, width_as)
        filename = urlIm
        
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        filename = FitsOutputPath + "{}_{}-{}_{:.0f}arcsec.fits.gz".format(radec_str, "2MASS",\
                                                                            params["fname"][:-8], width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename

    

# -------------------------------------------------------------------------- #
#                                   CFHT                                     #
# -------------------------------------------------------------------------- #
def cfht_dl(ra, dec, band, instrument="MegaPrime", optFiltersGen='both', getFullIm=False, width_as=20.,\
            FitsOutputPath="/data/cfht/", saveFITS=False):
    """
        Download CFHT fits image
        http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/doc/data/
        http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/doc/tap/
        
        Data access:
            http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/{ArchiveName}/{filename}
                CFHT filename = {ProductID}.fits.fz
        
        1/ Get the product ID
            http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/AdvancedSearch/tap/sync?LANG=ADQL&REQUEST=doQuery&USEMAQ=true&QUERY={query}&FORMAT=csv
                ex of ADQL query: ADQL tab on http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/search/#queryFormTab
                
        2/ Access the cutout
            http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:{ArchiveName}/{filename}&cutout=Circle+ICRS+{ra}+{dec}+{radius_deg}
            
       ------------------------------------------------------------------------------------------------------------------------
       
        instrument: CFHT intrument (MegaPrime, WIRCam, CFH12K MOSAIC)
        optFiltersGen: version of the MegaPrime filters (1stGen, 3rdGen or both)
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        width_as: size of the cutout if getFullIm=False
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk (to FitsOutputPath/BRICKID[0:3]/BRICKID/) 
    """
    
    print("   Get the url of the fits file ...")
    
    ## Filters ID
    # MegaPrime
    if instrument == "MegaPrime":
        filterIDs_1st = {"u": "'u.MP9301'", "g": "'g.MP9401'", "r": "'r.MP9601'", "i": "'i.MP9701'", "z": "'z.MP9801'"}
        filterIDs_3rd = {"u": "'u.MP9302'", "g": "'g.MP9402'", "r": "'r.MP9602'", "i": "'i.MP9703'", "z": "'z.MP9901'"}
                     
        if optFiltersGen == "1stGen":
            filterID = filterIDs_1st[band]   
        elif optFiltersGen == "3rdGen":
            filterID = filterIDs_3rd[band]                  
        elif optFiltersGen == 'both':
            filterID = filterIDs_1st[band] + "," + filterIDs_3rd[band]
            
    # WIRCam        
    elif instrument == "WIRCam":
        filterIDs = {"Y": "'Y.WC8002'", "J": "'J.WC8101'", "H": "'H.WC8201'", "K": "'Ks.WC8302'"}
        filterID = filterIDs[band]
        
    # CFH12k MOSAIC
    elif instrument == "CFH12K MOSAIC":
        filterID = "'" + band + "'"
    
    ## radec_str
    radec_str = radecStr(ra, dec, precision=1)

    ## 1/ Get the product ID (1arcsec search around input position)
    query = """
    SELECT
	 Plane.productID AS "Product ID",
	 Plane.time_exposure AS "Int. Time",
	 Observation.instrument_name AS "Instrument",
	 Plane.energy_bandpassName AS "Filter",
	 Observation.proposal_pi AS "P.I. Name",
	 Plane.position_sampleSize AS "Pixel Scale",
	 isDownloadable(Plane.publisherID) AS "DOWNLOADABLE"
    FROM caom2.Plane AS Plane 
	JOIN caom2.Observation AS Observation ON Plane.obsID = Observation.obsID 
    WHERE  ( INTERSECTS( CIRCLE('ICRS',
	{},
	 {},
	 2.777777777777778E-4),
	 Plane.position_bounds ) = 1 
    AND Observation.collection = 'CFHT'
    AND Observation.instrument_name = '{}'
    AND Plane.energy_bandpassName IN ({})
    AND Plane.calibrationLevel = '2' 
    AND Plane.dataProductType = 'image' 
    AND Observation.type = 'OBJECT' 
    AND  ( Plane.quality_flag IS NULL OR Plane.quality_flag != 'junk' ) )""".format(str(ra), str(dec), instrument, filterID)
    
    if sys.version_info.major == 2:
        urlQuery = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/AdvancedSearch/tap/sync?LANG=ADQL&REQUEST=doQuery&USEMAQ=true&" +\
                    urllib.urlencode({"QUERY":query, "FORMAT":"csv"}) 
    else:
        urlQuery = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/AdvancedSearch/tap/sync?LANG=ADQL&REQUEST=doQuery&USEMAQ=true&" +\
                    urllib.parse.urlencode({"QUERY":query, "FORMAT":"csv"})     

    t = ascii.read(urlQuery, header_start=0, delimiter=',')
    # Available to download
    try:
        t = t[t['"DOWNLOADABLE"'].mask == False]
    except:
        pass
    
    ## 2/ Cutout url 
    radius_deg = width_as / 2.0 / 3600.
    #if instrument == "CFH12K MOSAIC":
    #    radius_deg = radius_deg * 50.
        
    
    # Check if pixel at source position 
    if len(t) > 0:
        # Sort by decreasing exposure time
        t.sort('"Int. Time"')
        t.reverse()
        #t = t[0:50]   ## First 50 lines only
        for i in range(len(t)):
            productID = t['"Product ID"'][i]
            urlIm = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:CFHT/{}.fits.fz&cutout=Circle+ICRS+{:f}+{:f}+{:f}".format(productID, ra, dec, radius_deg)
            fname = "{}_CFHT-{}-{}_{:.0f}arcsec.fits.fz".format(radec_str, productID, band, width_as)
            FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
            try:
                res = urllib2.urlopen(urlIm)
                break
            except:
                urlIm = ""
                continue
                
    elif len(t) == 0 or urlIm == "":
        print("NO CFHT {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return ""
    
    # Full image 
    if getFullIm:
        # multi-extension file (one extension / CCD)
        urlIm = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/CFHT/{}.fits.fz".format(productID)
        fname = "CFHT-{}-{}.fits.fz".format(productID, band)
       
    filename = urlIm
        
    ## Save fits    
    if saveFITS:
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        filename = FitsOutputPath + fname
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename
    
    
    
# -------------------------------------------------------------------------- #
#                                  DECaLS                                    #
# -------------------------------------------------------------------------- #
def decals_dl(ra, dec, band, dr=5, width_as=20., getFullIm=False, brickid=None,\
              FitsOutputPath="/data/vault/des/DECALS/DR5/coadd/", saveFITS=False):
    """
        Download DECaLS fits cutouts
        
        A/ Full tile
         1/ Get brickid 
            BRICKID from survey-bricks file
            data disk: /data/vault/des/DECALS/DR5/survey-bricks.fits.gz
            url: http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr5/survey-bricks.fits.gz
            
         2/ coadd url
            ex: http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr5/coadd/010/0101m002/legacysurvey-0101m002-image-g.fits.fz
    
            prefix url: http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr5/coadd/
            file      : prefix url + BRICKID[0:3]/BRICKID/FILENAME
                filename: legacysurvey-BRICKID-image-BAND.fits.fz

        B/ Cutouts
           ex: http://legacysurvey.org/viewer/fits-cutout/?ra=190.1086&dec=1.2005&layer=decals-dr3&pixscale=0.27&bands=grz&size=512
         
         ------------------------------------------------------------------------------------------------------------------------
        dr: DECalS data release, default is 5
        width_as: width of the cutout in arcsec (used only if getFullIm is False)
        getFullIm: get the full DECaLs tile if True (default is False)
        brickid: brickid of the DECaLS observation, i.e. in the form "2054p095"
                If None, get it from survey-bricks table
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk (to FitsOutputPath/BRICKID[0:3]/BRICKID/)
          
    """
    
    print("   Get the url of the fits file ...")
    
    ## convert width to pixels
    pixelscale = 0.262
    width_pix = int(round(width_as / pixelscale, 0))
    
    ## A/ Full DECaLS tile
    if getFullIm:
        # 1/ Get brick and brickid
        if brickid is None:
            url_surveyBricks = "http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr{:d}/survey-bricks.fits.gz".format(dr)
                
            t_surveyBricks = Table.read(url_surveyBricks)
            ids = np.where((ra < t_surveyBricks["RA2"]) & (ra > t_surveyBricks["RA1"]) &\
                            (dec < t_surveyBricks["DEC2"]) & (dec > t_surveyBricks["DEC1"]))[0]
            if len(ids) > 0:
                brickid = t_surveyBricks["BRICKNAME"][ids][0]
                brick = brickid[0:3]
            else:
                print("NO DECaLS {}-BAND COVERAGE FOR THIS POSITION".format(band))
                return "", "", ""
                
        else:
            brick = brickid[0:3]
        
        # 2/ Coadd
        url = "http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr{:5}/coadd/".format(dr) +\
                     brick + "/" + brickid + "/" + "legacysurvey-{}-image-{}.fits.fz".format(str(brickid), band) 
        filename =  url  
        
    ## B/ Cutouts
    else:
        brickid = None
        url = "http://legacysurvey.org/viewer/fits-cutout/?ra={:f}&dec={:f}&layer=decals-dr{:d}&pixscale={}&bands={}&size={}".format(ra,\
                      dec, dr, pixelscale, band, width_pix)
        try:
            datas = fits.getdata(url, 0)
            filename = url
        except:
            print("NO DECaLS {}-BAND COVERAGE FOR THIS POSITION".format(band))
            return "", None
        
                       
    if saveFITS:
        if getFullIm:
            fname = "{}-legacysurvey-{}-image-{}.fits.fz".format("DECaLS", str(brickid), band)
            FitsOutputPath = FitsOutputPath + "{}/{}/".format(brick, brickid)
        else:
            radec_str = radecStr(ra, dec, precision=1)
            fname = "{}-{}-legacysurvey-image-{}_{:.0f}arcsec.fits".format("DECalS", radec_str, band, width_as)
            
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + fname
        os.system("wget -O {} '{}' ".format(filename, url))


    return filename, brickid
    


# -------------------------------------------------------------------------- #
#                                    DES                                     #
# -------------------------------------------------------------------------- #
def des_dl(ra, dec, band, dr=1, width_as=20., FitsOutputPath="/data/des/", saveFITS=False):
    """
        Download DES fits cutouts
        
        Example: https://datalab.noao.edu/desdr1/analysis/DwarfGalaxyDESDR1_20171101.html
          -> Chapter 3 -- Retrieve images
         
         urlIm = http://datalab.noao.edu/svc/cutout?col=des_dr1&siaRef=DES0145-5540_r2624p02_z.fits.fz&extn=1&POS=25.79265,-55.75297&SIZE=0.00555555555556,0.00555555555556
         
         ------------------------------------------------------------------------------------------------------------------------
        dr: DES data release, default is 1
        width_as: width of the cutout in arcsec 
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk (to FitsOutputPath/BRICKID[0:3]/BRICKID/)
          
    """
    from pyvo.dal import sia
    from astropy.io import votable
    import warnings
    warnings.filterwarnings('ignore', category=Warning, append=True)
    
    print("   Get the url of the fits file ...")
    
    

    # Input parameters
    width_deg = width_as / 3600.
    
    ## Get the deepest stacked images
    DEF_ACCESS_URL = "https://datalab.noao.edu/sia/des_dr{:d}".format(dr)
    svc = sia.SIAService(DEF_ACCESS_URL)
    imgTable = svc.search((ra,dec), (width_deg/np.cos(dec*np.pi/180), width_deg), verbosity=2).votable.to_table()
    
    # Select band and stack images
    sel = (imgTable["obs_bandpass"].astype(str)==band) &\
          ((imgTable["proctype"].astype(str)=="Stack") & (imgTable["prodtype"].astype(str)=="image")) 
    t = imgTable[sel]
    
    if len(t) > 0:
        # Deepest image (i.e. longest exposure time)
        row = t[np.argmax(t["exptime"].data.data.astype("float"))] 
        urlIm = row["access_url"].decode() # get the download URL
        pixelscale = row["im_scale"].data[0]
        fname = row["obs_publisher_did"].split("/")[-1][:-3]
        filename = urlIm

    else:
        print("NO DES-DR{:d} {}-band COVERAGE FOR THIS POSITION".format(dr, band))
        return "", 0.27
        
    ## Save the fits image
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + "{}-{}_{:.0f}arcsec.fits".format(radec_str, fname[:-5], width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))    
        
    return filename, pixelscale


# -------------------------------------------------------------------------- #
#                                  FLS KPNO                                  #
# -------------------------------------------------------------------------- #
def flsKPNO_dl(ra, dec, width_as=20., FitsOutputPath="/data/fls_kpno/", saveFITS=False):
    """
        Download FLS KPNO fits cutouts
        
        xFLS & ELAIS-N1 fields
         
         ------------------------------------------------------------------------------------------------------------------------
        width_as: width of the cutout in arcsec 
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk (to FitsOutputPath/BRICKID[0:3]/BRICKID/)
          
    """
    print("   Get the url of the fits file ...")
    
    ## Input params
    if (ra > 240) and (ra < 245):
        # ELAIS-N1 field
        mission = "FLS_ELAISN1_R"
        tbl = "cutouttbl1=FLS_ELAISN1_R"

    elif (ra > 255) and (ra < 265):
        # xFLS field
        mission = "FLS"
        tbl = "cutouttbl5=kpno_R&ntable_cutouts=5"
                
    else:
        print("NO FLS-KPNO R-band COVERAGE FOR THIS POSITION")
        return ""
        
    ## Get the cutouts
    url_search = "https://irsa.ipac.caltech.edu/cgi-bin/Cutouts/nph-cutouts?mission={}&units=arcsecg&locstr={:f}+{:f}&sizeX={:.2f}&min_size=1&max_size=180&{}&mode=PI".format(mission, ra, dec, width_as, tbl)
    
    root = ET.parse(urllib2.urlopen(url_search)).getroot()
    
    allfits = []
    for txt in root.findall("images/cutouts/fits"):
        link = txt.text.strip()
        if "fullfield" not in link:
            allfits.append(link)
            
    if len(allfits) > 0:
        urlIm = allfits[0]
        filename = urlIm
    else: 
        print("NO FLS-KPNO R-band COVERAGE FOR THIS POSITION")
        return ""  
    
    ## Save the fits image
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        
        fname = filename.split("/")[-1].split("_",3)[-1]    
        filename = FitsOutputPath + "{}-{}_{}_{:.0f}arcsec.fits".format(radec_str, "FLS-KPNO", fname[:-5], width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename
        
        
        
# -------------------------------------------------------------------------- #
#                     Hawaii Hubble Deep Field North                         #
# -------------------------------------------------------------------------- #
def hdfn_dl(ra, dec, band, width_as=20., FitsOutputPath="/data/hdfn/", getFullIm=False, saveFITS=False):   
    """
        A/ Download HDFN fits images (full image only, 1fits/band)
        http://www.astro.caltech.edu/~capak/hdf/
        
        b/ Cutouts using astroquery.skyview
           
       ------------------------------------------------------------------------------------------------------------------------
        width_as: width of the cutout in arcsec 
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        getFullIm: get full fits images if True
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
    """
    from astroquery.skyview import SkyView
     
    print("   Get the url of the fits file ...")
    
    # Input
    if band == "V":
        band = "V0201"
    
    if getFullIm:
        band = band.upper()
        base_url = "ftp://ftp.astro.caltech.edu/pub/capak/images/"
        fname = "HDF.{}.fits.gz".format(band)
        urlIm = base_url + fname 
        
    else:
        width_pix = int(round(width_as / 0.3, 0))
        c = SkyCoord(ra*u.deg, dec*u.deg)
        urlIm = SkyView.get_image_list(position=c, survey=['Hawaii HDF {}'.format(band)], pixels=str(width_pix))
        
        if len(urlIm) > 0:
            urlIm = urlIm[0]
        else:
            print("NO HDFN {}-BAND COVERAGE FOR THIS POSITION".format(band))
            return ""
        
        radec_str = radecStr(ra, dec, precision=1)
        fname = "{}-{}_{}_{:.0f}arcsec.fits".format(radec_str, "HDFN", band, width_as)

    filename = urlIm
    
    ## Save the fits image
    if saveFITS:
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + fname
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename
    
    
# -------------------------------------------------------------------------- #
#                                    HSC                                     #
# -------------------------------------------------------------------------- #
def hsc_dl(ra, dec, band, width_as=20., login=[], FitsOutputPath="/data/hsc/", saveFITS=False):   
    """
        A/ Download HSC fits images 
        https://hsc-release.mtk.nao.ac.jp/das_quarry/#cwh
        
        ACCOUNT REQUIRED TO ACCESS PUBLIC DATA
        
        ex: urlIm = https://hsc-release.mtk.nao.ac.jp/das_quarry/cgi-bin/quarryImage?ra=214.72951549960177&dec=52.120033840566514&sw=10arcsec&sh=10arcsec&type=coadd&image=on&filter=HSC-I&tract=&rerun=any
       ------------------------------------------------------------------------------------------------------------------------
       
        width_as: width of the cutout in arcsec 
        login: array containing the username and password to access HSC public database
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
    """
    
    import requests, re
      
    print("   Get the url of the fits file ...")
    
    # Filters
    filters = {"g": "HSC-G", "r": "HSC-R", "i": "HSC-I", "z": "HSC-Z", "y": "HSC-Y", "N816": "NB0816", "N921": "NB0921"}
    filter = filters[band]
    
    # login
    username, password = login[0], login[1]
    
    ## Cutouts url
    urlIm = "https://hsc-release.mtk.nao.ac.jp/das_quarry/cgi-bin/quarryImage?ra={:f}&dec={:f}&sw={:.0f}arcsec&sh={:.0f}arcsec&type=coadd&image=on&filter={}&tract=&rerun=any".format(ra, dec, width_as/2, width_as/2, filter)    
    filename = urlIm
    
    # url exists & fname
    r = requests.get(urlIm, auth=(username, password))
    
    if r.status_code == 404:
        print("NO HSC {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return ""  
        
    else:    
        cd = r.headers.get('content-disposition')
        fname = re.findall('filename=(.+)', cd)[0]
        fname = fname[8:-1]
        
    ## HTTP AUthentification
    p = urllib2.HTTPPasswordMgrWithDefaultRealm()
    p.add_password(None, filename, username, password)
    handler = urllib2.HTTPBasicAuthHandler(p)
    opener = urllib2.build_opener(handler)
    urllib2.install_opener(opener)
    
    ## Save the fits image
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + "{}_{}_{:.0f}arcsec.fits".format(radec_str, fname[:-5], width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename
    
        
# -------------------------------------------------------------------------- #
#                             HST/ACS - HST/WFC3                             #
# -------------------------------------------------------------------------- # 
def hst_dl(ra, dec, band, instrument="ACS", width_as=20.,\
           FitsOutputPath="/data/hst/", saveFITS=False):   
    """
        Hubble Legacy Archive Search
        http://hla.stsci.edu/hla_help.html -> VO servive
        http://hla.stsci.edu/fitscutcgi_interface.html
        
        1/ get filename
        http://hla.stsci.edu/cgi-bin/hlaSIAP.cgi?POS=189.151338,62.092132&SIZE=0.001&inst=ACS&Spectral_Elt=f435w
        instr: WFC3, ACS / Spectral_Elt: filter
        
        2/ Cutout access
        http://hla.stsci.edu/cgi-bin/fitscut.cgi?red=hlsp_goods_hst_acs-wfc_north-mosaic_f435w_v2.0_drz&RA=189.151338&Dec=62.092132&size=400&format=fits
        
        ------------------------------------------------------------------------------------------------------------------------
        
        band: HST filter
        instrument: HST instrument (ACS or WFC3)
        width_as: size of the cutout (in arcsec)
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)

    """
    from astropy.io.votable import parse_single_table
    from astropy.nddata import Cutout2D
    import warnings
    warnings.filterwarnings('ignore', category=Warning, append=True)
    
    print("   Get the url of the fits file ...")
    
    # Input parameters
    pixelscales = {"ACS": 0.05, "WFC3": 0.09}
    pixelscale = pixelscales[instrument]
    width_pix = int(round(width_as / pixelscale, 0))
    
    
    ## 1/ Get the filename
    urlSearch = "http://hla.stsci.edu/cgi-bin/hlaSIAP.cgi?POS={:f},{:f}&SIZE=0.001&imagetype=hlsp,mosaic,combined&inst={}&Spectral_Elt={}".format(ra, dec, instrument, band)
    t = parse_single_table(urlSearch).to_table()
    
    if len(t) > 0:
        # Data product rank: 5 (HLSP), 3 (mosaic), 2 (combined)
        t.sort(["Level"])
        t.reverse()
        
        # filenames
        fnames = t['filename']
       
    else:
        print("NO HST/{} {} COVERAGE FOR THIS POSITION".format(instrument, band))
        return ""
        
    ## 2/ Cutout access
    partialIm_lens, partialIm_fnames = [], []
    filename = ""
    
    for fname in fnames:
        fname = fname.decode('utf-8')   #Python3: convert bytes to str
        urlIm = "http://hla.stsci.edu/cgi-bin/fitscut.cgi?red={}&RA={:f}&Dec={:f}&size={:d}&format=fits".format(fname, ra, dec, width_pix)
        data = fits.getdata(urlIm)
        # Check if source fully covered by the image
        if np.isfinite(np.mean(data)):
            filename = urlIm
            break
        else:
            small_cutout = Cutout2D(data, (width_pix/2, width_pix/2), (0.01*width_pix, 0.01*width_pix))
            image = small_cutout.data
            if np.isfinite(np.nanmean(image)):
                partialIm_lens.append(len(data[np.isfinite(data[0])]))
                partialIm_fnames.append(fname)
                
    if filename == "" and len(partialIm_fnames)>0:
        filename = partialIm_fnames[np.argmax(partialIm_lens)]
            

    ## Save the fits image
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + "{}-{}_{}_{:.0f}arcsec.fits".format(radec_str, "HST", fname[:-5], width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename
    

    
# -------------------------------------------------------------------------- #
#                               Legacy Survey                                #
# -------------------------------------------------------------------------- #
def legacySurvey_dl(ra, dec, band, dr=6, width_as=20.,\
              FitsOutputPath="/data/LegacySurvey/", saveFITS=False):
    """
        Download Legacy Survey (BASS+MzLS) fits cutouts
        
        Cutouts
           ex: http://legacysurvey.org/viewer/fits-cutout?ra=154.7709&dec=46.4537&layer=mzls+bass-dr6&pixscale=0.27&bands=grz&size=512
         
         ------------------------------------------------------------------------------------------------------------------------
        dr: LegacySurvey data release, default is 6
        width_as: width of the cutout in arcsec (used only if getFullIm is False)
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk 
          
    """
    
    print("   Get the url of the fits file ...")
    
    ## convert width to pixels
    pixelscale = 0.262
    width_pix = int(round(width_as / pixelscale, 0))
    
    ## Cutouts
    url = "http://legacysurvey.org/viewer/fits-cutout/?ra={:f}&dec={:f}&layer=mzls+bass-dr{:d}&pixscale={}&bands={}&size={}".format(ra,\
                      dec, dr, pixelscale, band, width_pix)
                      
    try:
        datas = fits.getdata(url, 0)
        filename = url
    except:
        print("NO Legacy Survey {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return ""
                 
                 
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"     
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        fname = "{}-{}-legacysurvey-image-{}_{:.0f}arcsec.fits".format("LegacySurvey", radec_str, band, width_as)    
        filename = FitsOutputPath + fname
        os.system("wget -O {} '{}' ".format(filename, url))

    return filename
    
    
# -------------------------------------------------------------------------- #
#                                   NDWFS                                    #
# -------------------------------------------------------------------------- #
def ndwfs_dl(ra, dec, band, width_as=20., FitsOutputPath="/data/ndwfs/", saveFITS=False):
    """
        Download NDWFS fits cutouts
        
        http://r2.sdm.noao.edu/ndwfs/data-cutout.html
        filters = Bw, R, I, K
         
         urlIm = http://r2.sdm.noao.edu/ndwfs/cutout.php?ra=14:28:07&dec=34:55:40&rawidth=0.33&decwidth=0.33&filters=K
         
         ------------------------------------------------------------------------------------------------------------------------
        width_as: width of the cutout in arcsec 
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits tile file on disk (to FitsOutputPath/BRICKID[0:3]/BRICKID/)
          
    """
    print("   Get the url of the fits file ...") 
    
    # Input parameters
    width_am = width_as / 60.
    
    radec = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    ra_hms = "{:.0f}:{:.0f}:{:f}".format(radec.ra.hms.h, radec.ra.hms.m, radec.ra.hms.s)
    dec_dms = "{:.0f}:{:.0f}:{:f}".format(radec.dec.dms.d, radec.dec.dms.m, radec.dec.dms.s)
    
    ## Get the url of the cutout
    urlIm = "http://r2.sdm.noao.edu/ndwfs/cutout.php?ra={}&dec={}&rawidth={:.2f}&decwidth={:.2f}&filters={}".format(ra_hms, dec_dms,\
             width_am, width_am, band)
             
    try:
        h = fits.open(urlIm)
        h.close()
        filename = urlIm
    except:
        print("NO NDWFS {}-band COVERAGE FOR THIS POSITION".format(band))
        return ""
     
     
    ## Save the fits image
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + "{}-NDWFS_{}_{:.0f}arcsec.fits".format(radec_str, band, width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))    
        
    return filename


    
# -------------------------------------------------------------------------- #
#                               Pan-STARRS                                   #
# -------------------------------------------------------------------------- #    
def ps1_dl(ra, dec, band, radec_str, projcell=None, subcell=None, getFullIm=False, width_as=20.,\
           FitsOutputPath="/data/ps1ardata/ps1dr1_cutouts_cache/", saveFITS=False):
    """
        Download PanSTARRS fits file:
        https://outerspace.stsci.edu/display/PANSTARRS/PS1+Image+Cutout+Service
    
        1/ Full tile
            ex: http://ps1images.stsci.edu/FILENMAE
        2/ Cutout: size=int (in pixels) / imagename = output filename
            ex: http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=FILENAME&RA={}&DEC={}&size={}&format=fits&imagename={}
    
            FILENAME: /rings.v3.skycell/PROJCELL/SUBCELL/rings.v3.skycell.PROJCELL.SUBCELL.stk.BAND.unconv.fits
                = /rings.v3.skycell/PROJCELL/SUBCELL/SHORTNAME
            
            PROJCELL, SUBCELL from
            http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=" + str(ra) + "&dec=" + str(dec)
                -> return projcell subcell ra dec filter mjd type filename shortname
            (one line per band)
         
         ------------------------------------------------------------------------------------------------------------------------
    
        projcell, subcell: projection id and skycell id of the PS1 observation; length=4, 3 respectively (ex 1070, 079)
                If None, get it from http://ps1images.stsci.edu/cgi-bin/ps1filenames.py
        getFullIm: get the full fits image (if True) or a fits cutout of size width_as (if False; default)
        width_as: size of the cutout if getFullIm=False
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
          
    """
    
    print("   Get the url of the fits file ...")
    
    ## convert width to pixels
    pixelscale = 0.258
    width_pix = int(round(width_as / pixelscale, 0))
       
    ## Get skycellid and projectionid
    if projcell == None:
        url_filenames = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=" + str(ra) + "&dec=" + str(dec)
        datas = urllib2.urlopen(url_filenames).readlines()
        
        if len(datas) > 1:
            spilt_datas = datas[1].split(" ")
            projcell, subcell = spilt_datas[0], spilt_datas[1]
            projcell, subcell = projcell.rjust(4, "0"), subcell.rjust(3, "0")
        else:
            print("NO PanSTARRS {}-BAND COVERAGE FOR THIS POSITION".format(band))
            return "", None, None
        
    ## filename and shortname
    shortname = "rings.v3.skycell.{}.{}.stk.{}.unconv.fits".format(projcell, subcell, band)
    filename = "/rings.v3.skycell/{}/{}/{}".format(projcell, subcell, shortname)
        
    ## Get the fits image
    if getFullIm:
        # Get the full fits image
        url_im = "http://ps1images.stsci.edu/" + filename
    else:
        # Get a cutout of size width_as
        output_filename = radec_str + "_" + shortname
        url_im = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red={}&RA={}&DEC={}&size={}&format=fits&imagename={}".format(filename, str(ra), str(dec), str(width_pix), output_filename)

       
    ## Save the fits file    
    if saveFITS:
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        if getFullIm:
            filename = FitsOutputPath + "{}_{}_{}".format(radec_str, "PS1", shortname)
        else:
            filename = FitsOutputPath + "{}_{}-{}_{:.0f}arcsec.fits".format(radec_str, "PS1",\
                                                                        shortname[:-5], width_as)
        os.system("wget -O {} '{}'".format(filename, url_im[0]))
        url_im = filename
       
    return url_im, projcell, subcell
    
    
    
# -------------------------------------------------------------------------- #
#                                    SDSS                                    #
# -------------------------------------------------------------------------- #   
def sdss_dl(ra, dec, band, dr=12, objid=None, FitsOutputPath="/data/sdss/", saveFITS=False):   
    """
        Download SDSS fits images (entire plate only)
        http://skyserver.sdss.org/dr12/en/help/docs/api.aspx#imgcutout
        
        1/ Get RERUN, RUN, CAMCOL, and FIELD
           http://skyserver.sdss.org/dr{}/SkyserverWS/SearchTools/SqlSearch
        
        2/ Get the url of the fits image
           http://data.sdss3.org/sas/dr{dr}/boss/photoObj/frames/{rerun}/{run}/{camcol}/frame-{band}-{run:06d}-{camcol}-{field:04d}.fits.bz2
           
       ------------------------------------------------------------------------------------------------------------------------
    
        dr: SDSS data release (default is 12)
        objid: SDSS objid of the source (used if not None to query RERUN, RUN, CAMCOL, and FIELD; otherwise use position)
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
    """
    
    print("   Get the url of the fits file ...")
    
    ## 1/ Get RERUN, RUN, CAMCOL, and FIELD
    sqlQuery_prefix = "SELECT distinct p.run, p.rerun, p.camcol, p.field FROM photoobjall AS p"
                    
    if objid is None:
        sqlQuery = sqlQuery_prefix + " JOIN dbo.fGetNearestObjEq({:f}, {:f}, 1.0) AS f ON p.objid=f.objid".format(ra, dec)
    else:
        sqlQuery = "SELECT run,rerun,camcol,field FROM photoobjall WHERE objid={}".format(objid)
    
    if sys.version_info.major == 2:
        urlQuery = "http://skyserver.sdss.org/dr{}/SkyserverWS/SearchTools/SqlSearch?{}".format(dr,\
                                   urllib.urlencode({'cmd':sqlQuery,'format':'csv'}))
    else:
        urlQuery = "http://skyserver.sdss.org/dr{}/SkyserverWS/SearchTools/SqlSearch?{}".format(dr,\
                                   urllib.parse.urlencode({'cmd':sqlQuery,'format':'csv'}))
    t = ascii.read(urlQuery)
    
    if len(t) > 0:
        t = t[0]
        params = {'run': t["run"], 'rerun': t["rerun"], 'camcol': t["camcol"], 'field': t["field"]}
        
        ## 2/ Get the url of the fits image
        base_url = "http://data.sdss3.org/sas/dr{:d}/boss/photoObj/frames/{rerun}/{run}/{camcol}/".format(dr, band, **params)
        fname = "frame-{}-{run:06d}-{camcol}-{field:04d}.fits.bz2".format(band, **params)
        urlIm = base_url + fname
        params['base_url'] = base_url
        filename = urlIm
        
    else:
        print("NO SDSS {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return "", {}
            
    ## Save the fits image
    if saveFITS:
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + "{}-{}".format("SDSS", fname)
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename, params
    

 
# -------------------------------------------------------------------------- #
#                                  Spitzer                                   #
# -------------------------------------------------------------------------- #   
def spitzer_dl(ra, dec, band, width_as=20., dataset="SEIP", FitsOutputPath="/data/spitzer/", saveFITS=False):
    """
        Download Spitzer / IRAC fits images
        http://sha.ipac.caltech.edu/applications/Spitzer/SHA/help/doc/api.html
    
        1/ Get cutouts search results (XML output -> mode=PI)
            https://irsa.ipac.caltech.edu/cgi-bin/Cutouts/nph-cutouts?mission=SEIP&units=arcsecg&locstr={ra}{+/-dec}&sizeX={width_as}&min_size=18&max_size=1800&cutouttbl1=science&mode=PI
    
        2/ Get url for requested filter
            ex: https://irsa.ipac.caltech.edu:443/workspace/TMP_i8opJ4_24815/Cutouts/189.407970_62.2751050331.v0001/results/science/0001_189.40797000_62.27511000_70101860.70101860-0.IRAC.1.mosaic.fits
            
        NB: IRAC-EGS
        IRSA webpage -> only deep IRAc1, IRAC2 images/cutouts available from the s-CANDELS survey
        All bands images (full, 1 fits / band): https://www.cfa.harvard.edu/irac/egs/
        -------------------------------------------------------------------------------------------------------------    
                        
        width_as: size of the cutout in arcsec
        dataset: Spitzer dataset -> SEIP (default; includes SDWFS, SWIRE, SERVS), SDWFS, SERVS, SHELA, SIMES, SSDF, SpIES, SpUDS, SWIRE, GOODS 
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
    """
    
    print("   Get the url of the fits file ...")
    
    bandsNum_dict = {"I1": 1, "I2": 2, "I3": 3, "I4": 4}
    bandsNum = bandsNum_dict[band]
    width_deg = width_as / 3600.
    
    ## IRSA cutouts service
    IRSAcutouts_dict = {"SEIP": "Cutouts", "SDWFS": "Atlas", "SERVS": "Cutouts", "SHELA": "Cutouts", "SSDF": "Cutouts",\
                        "SpIES": "Atlas", "SpUDS": "Atlas", "SWIRE": "Cutouts", "GOODS": "Atlas", "FLS": "Cutouts"}
    IRSAcutouts = IRSAcutouts_dict[dataset]
    IRSAcutouttbl1_dict = {"SEIP": "science", "SDWFS": None, "SERVS": "irac", "SHELA": "cutouts", "SSDF": "irac",\
                           "SpIES": "cutouts", "SpUDS": "irac_sci", "SWIRE": "SWIRE_Spitzer_IRAC", "GOODS": "Spitzer_IRAC_DR3",\
                           "FLS": "irac_mosaics"}
    IRSAcutouttbl1 = IRSAcutouttbl1_dict[dataset]
    IRSAcutMin_dict = {"SEIP": "18", "SDWFS": None, "SERVS": "1", "SHELA": "1", "SSDF": "1", "SpIES": "1",\
                       "SpUDS": "756", "SWIRE": "1", "GOODS": "1", "FLS":"1"}
    IRAScutMin = IRSAcutMin_dict[dataset]
    IRSAcutMax_dict = {"SEIP": "1800", "SDWFS": None, "SERVS": "600", "SHELA": "601", "SSDF": "600", "SpIES": "601",\
                       "SpUDS": "7200", "SWIRE": "600", "GOODS": "600", "FLS":"180"}
    IRAScutMax = IRSAcutMax_dict[dataset]
    
    ## Filenames keyword
    fnameImType_dict = {"SEIP": ".mosaic.", "SDWFS": "combined_epochs", "SERVS": "mosaic", "SHELA": "s_wavg", "SSDF": "mosaic.",\
                        "SpIES": "Combined", "SpUDS": "mosaic", "SWIRE": "mosaic", "GOODS": "sci", "FLS": "mosaic"}
    fnameImType = fnameImType_dict[dataset]
    fnameBand_dict = {"SEIP": "IRAC.{:d}".format(bandsNum), "SDWFS": band, "SERVS": "irac.b{:d}".format(bandsNum),\
                      "SHELA": "ch{:d}".format(bandsNum), "SSDF": band, "SpIES": "ch{:d}".format(bandsNum),\
                      "SpUDS": "ch{:d}".format(bandsNum), "SWIRE": band, "GOODS": "irac_{:d}".format(bandsNum),\
                      "FLS": "chan{:d}".format(bandsNum)}
    fnameBand = fnameBand_dict[dataset]


    ## 1/ Get cutouts search results
    if IRSAcutouts == "Cutouts":
        xmlTags = "images/cutouts/fits"
        url_search = "https://irsa.ipac.caltech.edu/cgi-bin/Cutouts/nph-cutouts?mission={}&units=arcsecg&locstr={:f}+{:f}&sizeX={:.2f}&min_size={}&max_size={}&cutouttbl1={}&mode=PI".format(dataset, ra, dec, width_as, IRAScutMin, IRAScutMax, IRSAcutouttbl1)
 
    elif IRSAcutouts == "Atlas":
        url_search = "https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas?mission={}&locstr={:f}+{:f}&regSize={:f}&covers=on&mode=PI".format(dataset, ra, dec, width_deg)
        if dataset == "GOODS":
            url_search = "https://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas?mission=SGOODS&locstr={:f}+{:f}&regSize={:f}&covers=on&radius=0.05&radunits=deg&searchregion=on&mode=PI".format(ra, dec, width_deg)
        xmlTags = "images/metadata"
    
    root = ET.parse(urllib2.urlopen(url_search)).getroot()
    
    allfits = []
    for txt in root.findall(xmlTags):
        link = txt.text.strip()
        
        # Output link = tbl table or fits images
        if link.endswith(".tbl"):
            # Get fits url from output table (science data type only)
            t = Table.read(link, format='ipac')
            t = t[t["file_type"]=="science"]
            for res in range(len(t)):
                fname = t["fname"][res]
                # Cutout image
                urlIm = "https://irsa.ipac.caltech.edu/cgi-bin/Subimage/nph-subimage?origfile=/irsadata/SPITZER/{}/{}&ra={:f}&dec={:f}&xsize={:f}".format(dataset, fname, ra, dec, width_deg)
                if dataset == "SERVS":
                    # Only Full image available
                    urlIm = "https://irsa.ipac.caltech.edu/data/SPITZER/SERVS/{}".format(fname)
                
                ## 2/ Only mosaic image for requested filter
                if fnameBand in fname and fnameImType in fname:
                    allfits.append(urlIm)
            
        elif link.endswith(".fits"):
            ## 2/ Only mosaic image for requested filter
            if fnameBand in link and fnameImType in link:
                allfits.append(link)

    ## Get the fits mosaic image
    if len(allfits) > 0:
        urlIm = allfits[0]
        filename = urlIm
        
    else:
        print("NO Spitzer {} {}-BAND COVERAGE FOR THIS POSITION".format(dataset, band))
        return ""
    
    ## Save the fits image
    if saveFITS:
        if IRSAcutouts == "Cutouts":
            fname = filename.split("/")[-1].split("_",3)[-1]
        elif IRSAcutouts == "Atlas":
            fname = filename.split("/")[-1].split("&")[0]
            if dataset == "SpIES":
                fname = band + "_Combined_" + fname
            
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
            
        filename = FitsOutputPath + "{}_{}-{}_{}_{:.0f}arcsec.fits".format(radec_str, "Spitzer", dataset,\
                                                                            fname[:-5], width_as)
        os.system("wget -O {} '{}' ".format(filename, urlIm))
        
    return filename
  
  
  
  
# -------------------------------------------------------------------------- #
#                                    UHS                                     #
# -------------------------------------------------------------------------- #
def uhs_dl(ra, dec, band, wsaLogin=[], database="UHSDR1", width_as=20.,\
              FitsOutputPath="/data/uhs/", saveFITS=False):
    """
        Download UHS fits image
        
        1/ GetImage cut-out results
        ex: http://wsa.roe.ac.uk:8080/wsa/GetImage?database=UHSDR1&ra=104.7666&dec=16.4857&sys=J&name=&filterID=all&xsize=1&ysize=1&obsType=object&frameType=stack&mfid=&fsid=
        
        2/ Show GetImage cut-out:
        http://wsa.roe.ac.uk/cgi-bin/getImage.cgi?file=/disk40/wsa/ingest/fits/20130304_v5/w20130304_00629_st.fit&mfid=6228936&extNo=3&lx=299&hx=345&ly=397&hy=442&rf=90&flip=1&uniq=4930_730_9_77797_1&xpos=24&ypos=23.5&band=J&ra=104.7666&dec=16.4857
        
        3/ Get cut-out
        http://wsa.roe.ac.uk/cgi-bin/getFImage.cgi?file=/disk40/wsa/ingest/fits/20130304_v5/w20130304_00629_st.fit&mfid=6228936&extNo=3&lx=299&hx=345&ly=397&hy=442&rf=90&flip=1&uniq=4930_730_9_77797_1&xpos=24&ypos=23.5&band=J&ra=104.7666&dec=16.4857
        
        filename ex: w20130304_00629_st.fit
        
        ------------------------------------------------------------------------------------------------------------------------
        
        wsaLogin: array of [username, password] to login (Not compulsory, data puplic since 01/08/2018)
        database: UHS database used "UHSDR{}" (default if DR1)
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
        """

    print("   Get the url of the fits file ...")
    
    ## Login to UHS (Not compulsory, data puplic since 01/08/2018)
    if len(wsaLogin) !=0:
        login_url = "http://wsa.roe.ac.uk:8080/wsa/DBLogin?community=+&user={}&passwd={}&community2=UHS".format(wsaLogin[0], wsaLogin[1])
        response=urllib2.urlopen(login_url)
        request=urllib2.Request(login_url)
        cj=cookielib.CookieJar()
        cj.extract_cookies(response, request)
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))

    ## 1/ GetImage cut-out results
    width_am = width_as / 60.
    getImage = "http://wsa.roe.ac.uk:8080/wsa/GetImage?database={}&ra={}&dec={}&sys=J&name=&filterID=all&xsize={}&ysize={}&obsType=object&frameType=stack&mfid=&fsid=".format(database, ra, dec, width_am, width_am)
    
    ## 2/ Show GetImage cut-out
    if len(wsaLogin) !=0:
        res = opener.open(getImage).read()
    else:
        res = urllib2.urlopen(getImage).read()

    links = re.findall('href="(http://.*?)"', res.decode('utf-8'))
        
    if len(links) > 0:
        # 3/ Get cut-out
        filename = str(links[0])
        filename = filename.replace("getImage", "getFImage", 1)
    else:
        print("NO UHS {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return ""
            
    ## Save the fits file
    if saveFITS:
        split_url = filename.split("&")[0].split("/")
        fname = split_url[-1]
        
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        filename = FitsOutputPath + "{}_UHS-{}_{}_{:.0f}arcsec.fits.gz".format(radec_str, fname[:-4], band,\
                                                                            width_as)
        os.system("wget -O {} '{}'".format(filename, url_im[0]))
            
    return filename




# -------------------------------------------------------------------------- #
#                                  UKIDSS                                    #
# -------------------------------------------------------------------------- #
def ukidss_dl(ra, dec, band, wsaLogin=[], programme="LAS", database="UKIDSSDR10PLUS", width_as=20.,\
              FitsOutputPath="/data/ukidss/", saveFITS=False):
    """
        Download UKIDSS fits image
    
        from http://surveys.roe.ac.uk/wsa/cgi-bin/getFImage.cgi? using astroquery.ukidss
    
        ex: http://wsa.roe.ac.uk/cgi-bin/getFImage.cgi?file=/disk16/wsa/ingest/fits/20070220_v1/w20070220_02378_sf_st.fit&mfid=1239379&extNo=3&lx=3278&hx=3378&ly=439&hy=539&rf=90&flip=1&uniq=2936_987_16_80033_1&xpos=51.1&ypos=51&band=J&ra=205.5335707&dec=9.477173
         
            filename ex: 2936_987_16_80033_1.fits.gz
            
        ------------------------------------------------------------------------------------------------------------------------
    
        wsaLogin: array of [username, password, community] to access non-public data
        programme: survey name, i.e DXS, GCS, GPS, LAS, UDS (default is LAS)
        database: UKIDSS database used "UKIDSSDR{}PLUS" (default if DR10)
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
    """
    
    from astroquery.ukidss import Ukidss
    
    print("   Get the url of the fits file ...")
    
    ## Login to access non-public data
    if len(wsaLogin) !=0:
        u_obj = Ukidss(username=wsaLogin[0], password=wsaLogin[1], community=wsaLogin[2],\
                       database=database, programme_id=programme)
    else:
        u_obj = Ukidss(database=database, programme_id=programme)
                       
    ## Get the fits image
    width_am = width_as / 60.
    pos = SkyCoord(ra, dec, unit="deg")
    url_im = u_obj.get_image_list(pos, waveband=band, frame_type="stack", image_width=width_am*u.arcmin)
    
    if len(url_im) > 0:
        filename = url_im[0]
    else:
        print("NO UKIDSS {}-BAND COVERAGE FOR THIS POSITION".format(band))
        return ""
    
    ## Save the fits file    
    if saveFITS:
        split_url = filename.split("&")
        for string in split_url:
            if string.startswith("uniq"):
                fname = string[5:]

        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        filename = FitsOutputPath + "{}_{}-{}_{}_{:.0f}arcsec.fits.gz".format(radec_str, "UKIDSS", fname, band, width_as)
        os.system("wget -O {} '{}'".format(filename, url_im[0]))
        
    return filename



# -------------------------------------------------------------------------- #
#                                   VISTA                                    #
# -------------------------------------------------------------------------- #
def vista_dl(ra, dec, band, wsaLogin=[], survey="VHS", database="VHSDR5", width_as=20.,\
              FitsOutputPath="/data/vista/", saveFITS=False):
    """
        Download VISTA fits image (VHS, VIDEO, VIKING surveys)
        
        1/ GetImage cut-out results
        ex: http://horus.roe.ac.uk:8080/vdfs/GetImage?archive=VSA&database=VHSDR5&ra=152.0&dec=-20.&sys=J&filterID=all&xsize=0.33&ysize=0.33&obsType=object&frameType=tilestack&mfid=&fsid= 
           filterID: Z=1, Y=2, J=3, H=4, Ks=5 
           
        2/ Show GetImage cut-out:
        'http://horus.roe.ac.uk/wsa/cgi-bin/getImage.cgi?file=/disk35/vsa/ingest/fits/20130130_v1.3/v20130130_00563_st_tl.fit&mfid=2696055&extNo=1&lx=11236&hx=11295&ly=10956&hy=11015&rf=180&flip=1&uniq=5216_125_11_1727_1&xpos=30.3&ypos=30&band=J&ra=152&dec=-20
        
        3/ Get cut-out
        http://horus.roe.ac.uk/wsa/cgi-bin/getFImage.cgi?file=%2Fdisk35%2Fvsa%2Fingest%2Ffits%2F20130130_v1.3%2Fv20130130_00563_st_tl.fit;mfid=2696055;extNo=1;lx=11236;hx=11295;ly=10956;hy=11015;rf=180;flip=1;uniq=4916_666_39_1726_1;xpos=30.3;ypos=30;band=J;ra=152;dec=-20
        
        filename ex: v20130130_00563_st_tl.fit
        
        ------------------------------------------------------------------------------------------------------------------------
        
        wsaLogin: array of [username, password, community] to access non-public data
        survey: VISTA survey (VHS, VIKING, VIDEO)
        database: VISTA database used = survey + DataRealease (survey + data for proprietary data)
            Last public Release (up to June 2018): VHSDR5, VIKINGDR4, VIDEODR5
            full list: http://horus.roe.ac.uk:8080/vdfs/VgetImage_form.jsp
        FitsOutputPath: path to save the output fits file (if saveFITS=True)
        saveFITS: save fits file on disk (to FitsOutputPath/radec_str[0:4]/)
        """

    print("   Get the url of the fits file ...")
    
    ## VSA programme & filters ID
    filterIDs = {"Z":1, "Y":2, "J":3, "H":4, "K":5}
    filterID = filterIDs[band]
    
    ## Login to VISTA access non-public data
    if len(wsaLogin) !=0:
        login_url = "http://horus.roe.ac.uk:8080/vdfs/DBLogin?archive=VSA&community=&user={}&passwd={}&community2={}".format(wsaLogin[0],\
                             wsaLogin[1], wsaLogin[2])
        response=urllib2.urlopen(login_url)
        request=urllib2.Request(login_url)
        cj=cookielib.CookieJar()
        cj.extract_cookies(response, request)
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))

    # 1/ GetImage cut-out results
    width_am = width_as / 60.
    getImage = "http://horus.roe.ac.uk:8080/vdfs/GetImage?archive=VSA&database={}&ra={:f}&dec={:f}&sys=J&filterID={:d}&xsize={:.2f}&ysize={:.2f}&obsType=object&frameType=tilestack".format(database, ra, dec, filterID, width_am, width_am)

        
    # 2/ Show GetImage cut-out
    if len(wsaLogin) != 0:
        res = opener.open(getImage).read()
    else:
        res = urllib2.urlopen(getImage).read()

    links = re.findall('href="(http://.*?)"', res.decode('utf-8'))
            
    if len(links) > 0:
        # 3/ Get cut-out
        filename = links[0]
        filename = filename.replace("getImage", "getFImage", 1)
    else:
        print("NO VISTA-{} {}-BAND COVERAGE FOR THIS POSITION".format(survey, band))
        return ""
                
    ## Save the fits file
    if saveFITS:
        split_url = filename.split("&")[0].split("/")
        fname = split_url[-1]
            
        radec_str = radecStr(ra, dec, precision=1)
        FitsOutputPath = FitsOutputPath + radec_str[0:4] + "/"
        if not os.path.exists(FitsOutputPath):
            os.mkdir(FitsOutputPath)
        filename = FitsOutputPath + "{}_VISTA-{}_{}_{:.0f}arcsec.fits.gz".format(radec_str, survey, fname[:-4], band,\
                                                                              width_as)
        os.system("wget -O {} '{}'".format(filename, url_im[0]))
                
    return filename



             
