#!/usr/bin/python
import sys
import glob
import datetime as dt
from astropy.io import fits


def date2mjd(date):
  origin = dt.datetime(1858,11,17)
  mjd = (date-origin).days + (date-origin).seconds/86400.0
  return mjd

def mjd2date(mjd):
  origin = dt.datetime(1858,11,17)
  date = origin + dt.timedelta(mjd)
  return date

def mjd2jd(date):
    return date + 2400000.5

def jd2mjd(date):
    return date - 2400000.5

def get_timerange_in_ididata(datafits):
    """Returns the time range of the data recorded in this particular FITS data.

    Input
    -----
    datafits : astropy.io.fits.fitsrec.FITS_rec
        Data from a FITS-IDI file (obtained from e.g. astropy.io.fits.open(idifile)[#].data)

    Output
    ------
    starttime : datetime
        Starting time of the first scan in the file (in UT).
    endtime : datetime
        Ending time of the last scan in the file (in UT).
    """
    starttime = mjd2date(jd2mjd(datafits['DATE'][0] + datafits['TIME'][0]))
    endtime = mjd2date(jd2mjd(datafits['DATE'][-1] + datafits['TIME'][-1]))
    return starttime, endtime




def find_idi_with_time(idi_files, datetime=None, aipstime=None, verbose=True):
    """Returns the name of the IDI file that contains data corresponding to the given time.

    Input
    -----
    idi_files : list
        Path to all the FITS-IDI files. e.g. ['exp_1_1.IDI1', 'exp_1_1.IDI2', ..]
    datetime : datetime (optional, but either datetime or time must be provided)
        Date and time (in UT) to search for in the FITS-IDI files.
    aipstime : list (optional, but either datetime or time must be provided)
        Time to search for in the FITS-IDI files, in a AIPS time format (list with D, H, M, S).
    verbose : bool (default: True)
        Print messages.

    Output
    ------
    idi_file : str
        Path to the FITS-IDI file that contains the given time.
        It returns None if none of the files contain the given time.
    """
    # idi_files = glob.glob(path_idi)
    for idi_file in idi_files:
        if verbose:
            print('Opening {} file'.format(idi_file))
        hdu = fits.open(idi_file)
        hdu_data = hdu['UV_DATA'].data
        inittime, endtime = get_timerange_in_ididata(hdu_data)
        if datetime is not None:
            the_time = datetime
        elif aipstime is not None:
            the_time = dt.datetime.combine(inittime.date() + dt.timedelta(days=aipstime[0]),
                                           dt.time(*aipstime[1:]))
        else:
            raise Exception('Either datetime or time must be provided')

        if inittime < the_time < endtime:
            if verbose:
                print('{} contains {} UT'.format(idi_file, the_time.strftime('%d-%m-%Y %H:%M')))

            return idi_file

    if verbose:
        print('None of the files contain the requested time.')

    return None


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Two parameters are required!')
        print('%proc IDI_files YYYY/mm/dd/HH:MM')
        print('- IDI_files: path to the FITS-IDI files (wildcards accepted e.g. exp_1_1.IDI*)')
        print('- YYYY/mm/dd/HH:MM: time in UT to search for in the files.')
        sys.exit(1)

    the_idi_files = sys.argv[1:-1]
    if (len(the_idi_files) == 1) and ('*' in the_idi_files[0] or '?' in the_idi_files[1]):
        the_idi_files = glob.glob(the_idi_files[0])
    the_time = dt.datetime.strptime(sys.argv[-1], '%Y/%m/%d/%H:%M')
    find_idi_with_time(the_idi_files, the_time, verbose=True)


