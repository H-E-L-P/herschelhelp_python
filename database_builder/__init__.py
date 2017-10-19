# -*- coding: utf-8 -*-
#
# Script to add the field information to the database
#

import os
import sys
import warnings

from glob import glob

import numpy as np
from scipy.interpolate import CubicSpline

from astropy.io import votable
from astropy.table import Table
from astropy import units as u
from pymoc import MOC

HERE = os.path.dirname(__file__)
sys.path.append(os.path.join(HERE, '../'))
from herschelhelp.database import DATABASE_FILE, Database, Field, Filter


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


def add_filters():

    def _f99_extinction(wave):
        """Return interpolated Fitzpatrick 99 galactic extinction curve

        Parameters
        ----------
        wave : astropy.Quantity
            Input wavelength array, in units of length.

        Returns
        -------
        numpy array of floats
            Wavelength dependent extinction value.

        Notes
        -----

        Cubic spline anchorpoints taken from Table 3 of "Correcting for the
        Effects of Interstellar Extinction" by Fitzpatrick, Edward L.
        http://adsabs.harvard.edu/abs/1999PASP..111...63F

        """
        anchors_x = [0., 0.377, 0.820, 1.667, 1.828, 2.141, 2.433, 3.704,
                     3.846]
        anchors_y = [0., 0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265,
                     6.591]

        f99 = CubicSpline(anchors_x, anchors_y)
        inv_wave = (1 / wave.to(u.micron)).value

        return f99(inv_wave)

    with Database(writable=True) as database:
        for filter_file in glob("{}/filters/*.xml".format(HERE)):

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', Warning)
                filter_vo = votable.parse_single_table(filter_file)

            def _get_field(field_id):
                try:
                    field = filter_vo.get_field_by_id(field_id)
                except KeyError:
                    return None
                if field.unit is not None:
                    return field.value * field.unit
                elif field.datatype == 'char':
                    return field.value.decode('utf-8')
                else:
                    return field.value

            filter_id = os.path.basename(filter_file)[:-4]
            description = _get_field("Description")
            band_name = _get_field("Band")
            facility = _get_field("Facility")
            instrument = _get_field("Instrument")

            # Notes added by HELP
            try:
                notes = filter_vo.get_field_by_id("AdditionalProcessing")\
                    .description.replace("\n", "")
            except KeyError:
                notes = None

            # Filter response curve
            wave_unit = filter_vo.get_field_by_id("Wavelength").unit
            wave = filter_vo.array['Wavelength'].data * wave_unit
            transmission = filter_vo.array['Transmission'].data
            response = np.array([wave.to(u.angstrom).value, transmission])

            # Attenuation in the filter / E(B-V)
            att_ebv = (np.trapz(transmission * _f99_extinction(wave), wave) /
                       np.trapz(transmission, wave))
            att_ebv = np.round(att_ebv, 3)

            # Computing the mean wavelength because not all filters may
            # provide it.
            mean_wavelength = (
                np.trapz(response[1] * response[0], response[0]) /
                np.trapz(response[1], response[0]))

            # Computing the min and max wavelength with are the first and last
            # wavelength with at least 1% of the maximal transmission (we took
            # the definition from the Spanish VO site).
            min_wave_idx, max_wave_idx = np.where(
                response[1] > .01 * np.max(response[1])
            )[0][[0, -1]]
            min_wavelength = response[0][min_wave_idx]
            max_wavelength = response[0][max_wave_idx]

            database.session.add(
                Filter(filter_id, description, band_name, facility, instrument,
                       mean_wavelength, min_wavelength, max_wavelength,
                       att_ebv, response, notes)
            )

        database.session.commit()


def build_base():

    try:
        os.remove(DATABASE_FILE)
    except FileNotFoundError:
        pass

    with Database(writable=True) as database:
        database.upgrade_base()

    add_fields()
    add_filters()


if __name__ == "__main__":
    build_base()
