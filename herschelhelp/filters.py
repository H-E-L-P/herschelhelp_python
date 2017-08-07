# -*- coding: utf-8 -*-

from astropy.table import Table, Column

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

