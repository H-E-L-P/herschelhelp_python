# -*- coding: utf-8 -*-
#
# Script to add the field information to the database
#

import os
import sys

from astropy.table import Table
from pymoc import MOC

HERE = os.path.dirname(__file__)
sys.path.append(os.path.join(HERE, '../'))
from herschelhelp.database import DATABASE_FILE, Database, Field


def add_fields():

    all_moc = MOC()
    with Database(writable=True) as database:
        for field in Table.read("{}/fields.txt".format(HERE), format="ascii"):

            name = field['name']
            short_id = field['short_id']

            footprint = MOC()
            footprint.read("{}/coverages/{}_MOC.fits".format(HERE, name))

            database.session.add(Field(name, short_id, footprint))
            all_moc += footprint

        database.session.add(Field('__ALL__', '', all_moc))
        database.session.commit()


def build_base():

    try:
        os.remove(DATABASE_FILE)
    except FileNotFoundError:
        pass

    with Database(writable=True) as database:
        database.upgrade_base()

    add_fields()


if __name__ == "__main__":
    build_base()
