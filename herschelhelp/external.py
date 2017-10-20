# -*- coding: utf-8 -*-
import logging
import os

import numpy as np

from .filters import get_filter_meta_table

FILTER_MEAN_LAMBDAS = {
    item['filter_id']: item['mean_wavelength'] for item in
    get_filter_meta_table()
}
FILTER_MIN_LAMBDAS = {
    item['filter_id']: item['min_wavelength'] for item in
    get_filter_meta_table()
}

logging.basicConfig(level=logging.getLevelName(
    os.getenv('LOG_LEVEL', 'WARNING')))
LOGGER = logging.getLogger(__name__)


def convert_table_for_cigale(catalogue, inplace=False,
                             remove_zerofluxes=False):
    """Convert a HELP formated catalogue to be used by CIGALE.

    This function converts a HELP formated catalogue to a CIGALE formated one:

    - The “help_id” column is renamed to “id”;
    - Keep only the source having a redshift;
    - Only the total flux columns are kept plus some interesting columns;
    - The flux and error columns are renamed to <filter> and <filter>_err;
    - The fluxes that are flagged not to be used in SED fitting as set to NaN;
    - For each band, for the sources at a redshift for which the Lyman limit
      interfere with the filter bandpass, the flux is set to NaN.
    - If asked to remove them, set to NaN flux and error for fluxes equals
      to 0.
    - The fluxes are converted to mJy.

    Parameters
    ----------
    catalogue: astropy.table.Table
        The catalogue to convert.
    inplace: boolean
        If inplace is set to True, the function will not make a copy of the
        input catalogue.  This will save some memory space at the expense of
        modifying the input catalogue.
    remove_zerofluxes: boolean
        If True, all the 0 fluxes and associated error are set to NaN.

    Returns
    -------
    astropy.table.Table
        The converted catalogue.

    """
    if not inplace:
        catalogue = catalogue.copy()

    try:
        catalogue['help_id'].name = 'id'
    except KeyError:
        raise KeyError("The catalogue must have a help_id column.")

    try:
        catalogue = catalogue[~np.isnan(catalogue['redshift'])]
    except KeyError:
        raise KeyError("The catalogue must have a redshift column.")

    columns = ['id']
    # Keep interesting columns
    for col in ['ra', 'dec', 'ebv', 'redshift', 'stellarity', 'flag_gaia']:
        if col in catalogue.colnames:
            columns.append(col)

    # Total flux bands in the catalogue
    bands = [col[2:] for col in catalogue.colnames if col.startswith('f_')
             and not col.startswith('f_ap')]
    # Sort by mean wavelength
    bands.sort(key=lambda x: FILTER_MEAN_LAMBDAS.get(x, 0))

    for band in bands:

        LOGGER.debug("Processing band %s.", band)

        # Warn if the band is not in HELP database
        if band not in FILTER_MEAN_LAMBDAS:
            LOGGER.warning("The band %s is not in HELP filter database.",
                           band)

        # Flag telling not to use the band is present
        try:
            band_flag = catalogue['flag_{}'.format(band)]
        except KeyError:
            band_flag = None

        # Total flux
        LOGGER.debug("Converting %s flux to mJy.", band)
        catalogue['f_{}'.format(band)].name = band
        catalogue[band] /= 1000.  # μJy to mJy
        if band_flag is not None:
            catalogue[band][band_flag] = np.nan
        columns.append(band)

        # Error
        LOGGER.debug("Converting %s flux error to mJy.", band)
        if 'ferr_{}'.format(band) in catalogue.colnames:
            catalogue['ferr_{}'.format(band)].name = '{}_err'.format(band)
            catalogue['{}_err'.format(band)] /= 1000.  # μJy to mJy
            if band_flag is not None:
                catalogue['{}_err'.format(band)][band_flag] = np.nan
            columns.append('{}_err'.format(band))
        else:
            LOGGER.warning("The catalogue is missing the ferr_%s column.",
                           band)

        # 0 fluxes
        if remove_zerofluxes:
            LOGGER.debug("Remove %s fluxes equal to 0.", band)
            mask = ~np.isnan(catalogue[band])
            mask[mask] &= (catalogue[band][mask] == 0)
            catalogue[band][mask] = np.nan
            if '{}_err'.format(band) in catalogue.colnames:
                catalogue['{}_err'.format(band)][mask] = np.nan

        # Set to NaN the fluxes in the filters that overlap or are below the
        # Lyman break limit at the source redshift.
        if band in FILTER_MEAN_LAMBDAS:
            LOGGER.debug("Removing %s fluxes affected by Lyman limit.", band)
            lyman_limit_at_z = 912 * (1 + catalogue['redshift'])
            below_ly = np.where(FILTER_MIN_LAMBDAS[band] < lyman_limit_at_z)[0]
            if len(below_ly) > 0:
                catalogue[band][below_ly] = np.nan
                if 'ferr_{}'.format(band) in catalogue.colnames:
                    catalogue['{}_err'.format(band)][below_ly] = np.nan
                LOGGER.info("For %s sources, the band %s should not be used "
                            "because it overlaps or is below the Lyman limit "
                            "at the redshift of these sources. These fluxes "
                            "were set to NaN.",
                            len(below_ly), band)

    # Keep only the column we want
    return catalogue[columns]
