# -*- coding: utf-8 -*-

import logging
import os

from astropy.io.registry import IORegistryError
from astropy.table import Column, Table

from .database import get_filters

logging.basicConfig(level=logging.getLevelName(os.getenv('LOG_LEVEL', 'INFO')))
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
                  'instrument', 'mean_wavelength', 'att_ebv', 'notes']
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

    # TODO: Find a way to deal with not readable ascii files
    try:
        catalogue = Table.read(analysed_table)
    except IORegistryError:
        catalogue = Table.read(analysed_table, format='ascii')

    catalogue_bands = [col[5:] for col in catalogue.colnames
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
