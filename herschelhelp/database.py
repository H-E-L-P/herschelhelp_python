# -*- coding: utf-8 -*-

import logging
import os
import pkg_resources

from sqlalchemy import Column, PickleType, String, create_engine
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
    return depends of the arguments.

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

