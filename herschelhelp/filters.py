# -*- coding: utf-8 -*-

import logging
import os

from astropy.io import fits
from astropy.io.registry import IORegistryError
from astropy.table import Column, Table

from .database import get_filters

logging.basicConfig(level=logging.getLevelName(
    os.getenv('LOG_LEVEL', 'WARNING')))
LOGGER = logging.getLogger(__name__)


def get_filter_meta_table():
    """Generate a table with meta information about HELP filters

    This function generates an astropy.table.Table containing the information
    about the filters used in HELP, except their transmission profile.

    Returns
    -------
    astropy.table.Table

    """

    # List of filter attributes to put in the table
    attributes = ['filter_id', 'description', 'band_name', 'facility',
                  'instrument', 'mean_wavelength', 'min_wavelength',
                  'max_wavelength', 'att_ebv', 'notes']
    all_filters = get_filters()

    table_columns = []
    for attribute in attributes:
        table_columns.append(Column(
            name=attribute,
            data=[getattr(_, attribute) for _ in all_filters]
        ))

    return Table(table_columns)


def export_to_cigale(directory):
    """Export the HELP filter database to CIGALE compatible files.

    This function export the filter transmission profiles to files (one per
    filter) that can be imported in CIGALE and be used in SED fitting.

    Parameters
    ----------
    directory: str
        Directory in which to save the filter files.

    """
    for filt in get_filters():
        with open("{}/{}.dat".format(directory, filt.filter_id), 'w') as out:
            out.write("# {}\n".format(filt.filter_id))
            out.write("# energy\n")
            out.write("# {}".format(filt.description))
            Table(filt.response.T).write(out, format="ascii.no_header")


def export_to_eazy(analysed_table):
    """Export the HELP filter database to be used with EAZY.

    This function export the filter transmission profiles to be used by EAZY
    for photometric redshift computation.  As EAZY needs to associate to each
    column in the catalogue the index of the filter in the database, we are
    exporting it for a given table that will be analysed.  This table must be
    in the format used within HELP.

    It create three files <analysed_table>.res containing the responses of the
    filters in analysed_table, <analysed_table>.translate associating the
    aperture flux and error columns to their corresponding filter, and
    <analysed_table>.info containing additional information about the filters.

    Parameters
    ----------
    analysed_table: str
        The table that will be processed with EAZY. The aperture fluxes must be
        labelled f_ap_<filter> and the associated errors ferr_ap_<filter>.

    """

    response_filename = "{}.res".format(os.path.splitext(analysed_table)[0])
    translate_filename = "{}.translate".format(
        os.path.splitext(analysed_table)[0])
    info_filename = "{}.info".format(os.path.splitext(analysed_table)[0])

    if os.path.exists(response_filename):
        raise IOError("{} file exists.".format(response_filename))
    if os.path.exists(translate_filename):
        raise IOError("{} file exists.".format(translate_filename))
    if os.path.exists(info_filename):
        raise IOError("{} file exists.".format(info_filename))

    # If the table is a FITS file, we use astropy.io.fits so that we can deal
    # with huge files without loading them in memory.
    # TODO: Find a way to deal with not readable ascii files
    if analysed_table.endswith(".fits") or analysed_table.endswith(".fit"):
        with fits.open(analysed_table) as hdulist:
            column_names = hdulist[1].columns.names
    else:
        try:
            catalogue = Table.read(analysed_table)
        except IORegistryError:
            catalogue = Table.read(analysed_table, format='ascii')
        column_names = catalogue.colnames

    # EAZY uses aperture fluxes.
    catalogue_bands = [col[5:] for col in column_names
                       if col.startswith('f_ap_')]

    with open(response_filename, 'w') as filter_responses, \
            open(translate_filename, 'w') as translate_table, \
            open(info_filename, 'w') as info_file:

        for idx_band, band in enumerate(catalogue_bands):

            filt = get_filters(band)

            if filt is None:

                # The band is not in HELP filter database
                LOGGER.error("Filter %s is not in the database.", band)

            else:

                filt_nb_points = len(filt.response[0])
                filter_responses.write(
                    "{:>8d} {}, {}, mean_lambda={}, att/ebv={}\n".format(
                        filt_nb_points,
                        band,
                        filt.description,
                        filt.mean_wavelength,
                        filt.att_ebv))
                for idx, row in enumerate(filt.response.T):
                    filter_responses.write(
                        "{:>8d} {:>10.2f} {:>12.8g}\n".format(
                            idx+1, row[0], row[1])
                    )

                translate_table.write("f_ap_{} F{}\n".format(
                    band, idx_band + 1))
                translate_table.write("ferr_ap_{} E{}\n".format(
                    band, idx_band + 1))

                info_file.write(
                    "{:>3d} {}, {}, mean_lambda={}, att/ebv={}\n".format(
                        idx_band + 1,
                        band,
                        filt.description,
                        filt.mean_wavelength,
                        filt.att_ebv))


def correct_galactic_extinction(catalogue, inplace=False):
    """Correct photometric catalogue for galactic extinction.

    This function takes a photometric catalogue in the HELP format and correct
    the fluxes and magnitudes for the galactic extinction given the E(B-V)
    value associated with each source.  The catalogue must have an ebv column
    and the fluxes and magnitudes columns are identified using the HELP column
    format (f_<filter>, ferr_filter, f_ap_<filter>, ...).

    Column associated with filters that are not in the HELP database will not
    be corrected and an error message will be logged with the name of the
    missing filters.

    Parameters
    ----------
    catalogue: astropy.table.Table
        The catalogue to be corrected.
    inplace: boolean
        If inplace is set to True, the function will not make a copy of the
        input catalogue.  This will save some memory space at the expense of
        modifying the input catalogue.

    Returns
    -------
    astropy.table.Table
        The corrected catalogue.

    """
    if not inplace:
        catalogue = catalogue.copy()

    try:
        ebv = catalogue['ebv']
    except KeyError:
        raise KeyError("The catalogue is missing the ebv column.")

    # Instead of logging an error message for each missing band and problematic
    # column, lets just log a summary.
    missing_band, prob_columns = set(), set()

    for column in catalogue.colnames:

        band = None
        if column.startswith('f_ap_') or column.startswith('m_ap_'):
            band = column[5:]
        elif column.startswith('ferr_ap_') or column.startswith('merr_ap_'):
            band = column[8:]
        elif column.startswith('f_') or column.startswith('m_'):
            band = column[2:]
        elif column.startswith('ferr_') or column.startswith('merr_'):
            band = column[5:]

        if band is None:
            LOGGER.debug("Column %s is not associated with any band",
                         column)
        else:
            filt = get_filters(band)

            if filt is None:
                missing_band.add(band)
                prob_columns.add(column)
            else:
                att = filt.att_ebv * ebv

                if column.startswith("f"):
                    catalogue[column] *= 10**(att/2.5)
                elif column.startswith("merr"):
                    # Error in magnitude is not affected by extinction.
                    pass
                else:
                    catalogue[column] -= att

    if len(missing_band) > 0:
        LOGGER.error("The filters are missing in the database: %s.\n"
                     "The are present in these columns: %s.",
                     ", ".join(sorted(missing_band)),
                     ", ".join(sorted(prob_columns)))

    return catalogue
