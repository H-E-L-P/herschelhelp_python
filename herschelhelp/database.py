# -*- coding: utf-8 -*-

import logging
import os

import numpy as np

import pkg_resources

from sqlalchemy import Column, Float, PickleType, String, Text, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

logging.basicConfig(level=logging.getLevelName(os.getenv('LOG_LEVEL', 'INFO')))
LOGGER = logging.getLogger(__name__)

DATABASE_FILE = pkg_resources.resource_filename(__name__, "data.db")
ENGINE = create_engine('sqlite:///' + DATABASE_FILE, echo=False)
BASE = declarative_base()
SESSION = sessionmaker(bind=ENGINE)


class Field(BASE):
    """Sky field in HELP"""

    __tablename__ = 'fields'

    name = Column(String, primary_key=True)
    short_id = Column(String)
    footprint = Column(PickleType)

    def __init__(self, name, short_id, footprint):
        """Create a field.

        Parameters
        ----------
        name : string
            Name of the field.
        short_id : string
            Three letter identifiers of the field.
        footprint : pymoc.MOC
            Multi-order coverage map representing the footprint of the field.

        """
        self.name = name
        self.short_id = short_id
        self.footprint = footprint


class Filter(BASE):
    """Photometric filters

    This class represents a photometric filter in HELP database. A filter has
    these attributes:

    - filter_id: A unique identifier.  It is used in naming the column in HELP
      catalogues.
    - description: A human readable short description of the filter.
    - band_name: A standard representation for the band (e.g. “g” for a g band
      filter).
    - facility: The name of the observatory or telescope.
    - instrument: The name of the instrument.
    - mean_wavelength: The mean wavelength in Angstrom.
    - min_wavelength: First wavelength (in Angstrom) with a transmission at
      least 1% of the maximum transmission.
    - max_wavelength: Last wavelength (in Angstrom) with a transmission at
      least 1% of the maximum transmission.
    - att_ebv: The ratio between the attenuation in the filter and the E(B-V)
      color excess. It is computed using the Fitzpatrick (1999) extinction
      curve assuming a flat emission spectrum.
    - response: The filter response curve as a two axis numpy array with the
      wavelength in Angstrom in the first axis and the energetic transmission
      in the second one.
    - notes: Notes from the HELP project.

    """

    __tablename__ = 'filters'

    filter_id = Column(String, primary_key=True)
    description = Column(String)
    band_name = Column(String)
    facility = Column(String)
    instrument = Column(String)
    mean_wavelength = Column(Float)
    min_wavelength = Column(Float)
    max_wavelength = Column(Float)
    att_ebv = Column(Float)
    response = Column(PickleType)
    notes = Column(Text)

    def __init__(self, filter_id, description, band_name, facility, instrument,
                 mean_wavelength, min_wavelength, max_wavelength, att_ebv,
                 response, notes):
        """Create a photometric filter.

        Parameters
        ----------
        filter_id: string
            Unique identfier.
        description: string
            Human readable short description of the filter.
        band_name: string
            Standard representation of the band.
        facility: string
            Observatory or telescope.
        instrument: string
            Instrument name.
        mean_wavelength: float
            Mean wavelength of the filter in Angstrom.
        min_wavelength: float
            Minimal wavelength with at least 1% of maximum transmission in
            Angstrom.
        max_wavelength: float
            Last wavelength with at least 1% of maximum transmission in
            Angstrom.
        att_ebv: float
            Ratio between the attenuation in the filter and the E(B-V) color
            excess.
        response: numpy array of float
            Response curve of the filter, response[0] is the wavelength in
            Angstrom and transmission_curve[1] is the energetic transmission.
        notes: string
            Notes about the filter is any.
        """
        self.filter_id = filter_id
        self.description = description
        self.band_name = band_name
        self.facility = facility
        self.instrument = instrument
        self.mean_wavelength = mean_wavelength
        self.min_wavelength = min_wavelength
        self.max_wavelength = max_wavelength
        self.att_ebv = att_ebv
        self.response = response
        self.notes = notes


class Database(object):
    """Object giving access to the module database.

    This object can be used in a context manager.
    """

    def __init__(self, writable=False):
        """Create a connection manager to the module database

        writable : boolean
            If true, the users will have write access to the database
            provided they can write the sqlite file. By default, it's
            false.
        """
        self.session = SESSION()
        self.is_writable = writable
        LOGGER.debug("Initialising the database.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def upgrade_base(self):
        """Upgrade the table schemas in the database."""
        if self.is_writable:
            BASE.metadata.create_all(ENGINE)
        else:
            raise Exception("The database is not writable.")

    def close(self):
        """Close the connection to the database."""
        self.session.close_all()
        LOGGER.debug("Closing the database.")


def get_field(*args):
    """Retrieves HELP fields information.

    This function retrieves information about HELP fields in the database. The
    return depends on the arguments.

    - If there is only one argument, the function returns the field
      corresponding to the name given in argument. If there is no field with
      that name, None is returned.
    - If there are several arguments, the function return the list of fields
      corresponding to the given names. If none of the given names correspond
      to actual fields, an empty list is returned.
    - If no argument is given, the function returns the list of all the HELP
      fields.

    Each field is returned as a Field object with these attributes:

    - name: the full name of the field;
    - short_id: its three letters identifier;
    - footprint: its footprint as a pymoc.MOC object.

    Parameters
    ----------
    *args : str
        Name or names (several arguments, not a list of names) of the fields.

    Returns
    -------
    Field or list of fields
        If exactly one argument is given returns a Field (or None) else returns
        a list of fields.

    """

    with Database() as d:
        if len(args) == 0:
            return d.session.query(Field).all()
        elif len(args) == 1:
            return d.session.query(Field).get(args[0])
        else:
            return d.session.query(Field).filter(
                Field.name.in_(args)).all()


def get_filters(*args):
    """Retrieves HELP photometric filters information.

    This function retrieves information about the photometric filters used in
    HELP.  The return depends on the arguments.

    - If there is only one argument, the function returns the filter
      corresponding to the filter ID given in argument. If there is no filter
      with that ID, None is returned.
    - If there are several arguments, the function return the list of filters
      corresponding to the given IDs. If none of the given IDs correspond to
      actual filters, an empty list is returned.
    - If no argument is given, the function returns the list of all the filters
      in the database.

    Parameters
    ----------
    *args : str
        ID or IDs (several arguments, not a list) of the filters.

    Returns
    -------
    Filter or list of Filters
        If exactly one argument is given returns a Filter (or None) else
        returns a list of filters.

    """

    with Database() as d:
        if len(args) == 0:
            return d.session.query(Filter).all()
        elif len(args) == 1:
            return d.session.query(Filter).get(args[0])
        else:
            return d.session.query(Filter).filter(
                Filter.filter_id.in_(args)).all()
