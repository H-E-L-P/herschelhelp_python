# -*- coding: utf-8 -*-
import logging
import os

import numpy as np

from .filters import get_filter_meta_table

FILTER_MEAN_LAMBDAS = {
    item['filter_id']: item['mean_wavelength'] for item in
    get_filter_meta_table()
}

logging.basicConfig(level=logging.getLevelName(os.getenv('LOG_LEVEL', 'INFO')))
LOGGER = logging.getLogger(__name__)


def convert_table_for_cigale(catalogue, inplace=False):
    """Convert a HELP formated catalogue to be used by CIGALE.

    This function converts a HELP formated catalogue to a CIGALE formated one:

    - The “help_id” column is renamed to “id”;
    - Only the total flux columns are kept;
    - The flux and error columns are renamed to <filter> and <filter>_err;
    - The fluxes that are flagged not to be used in SED fitting as set to Nan;
    - The fluxes are converted to mJy.

    Parameters
    ----------
    catalogue: astropy.table.Table
        The catalogue to convert.
    inplace: boolean
        If inplace is set to True, the function will not make a copy of the
        input catalogue.  This will save some memory space at the expense of
        modifying the input catalogue.

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
        catalogue['f_{}'.format(band)].name = band
        catalogue[band] /= 1000.  # μJy to mJy
        if band_flag is not None:
            catalogue[band][band_flag] = np.nan
        columns.append(band)

        # Error
        if 'ferr_{}'.format(band) in catalogue.colnames:
            catalogue['ferr_{}'.format(band)].name = '{}_err'.format(band)
            catalogue['{}_err'.format(band)] /= 1000.  # μJy to mJy
            if band_flag is not None:
                catalogue['{}_err'.format(band)][band_flag] = np.nan
            columns.append('{}_err'.format(band))
        else:
            LOGGER.warning("The catalogue is missing the ferr_%s column.",
                           band)

    # Keep only the column we want
    return catalogue[columns]
