# -*- coding: utf-8 -*-

from astropy.table import Column, Table

from .database import get_filters


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
    """Export the HELP filter database to CIGALE compatible files

    This function export the filter transmission profiles to filtes (one per
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
