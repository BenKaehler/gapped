from gzip import GzipFile
from os.path import realpath, abspath, dirname, join
from inspect import getfile, currentframe

from cogent import Alignment

_data_dir = join(realpath(abspath(dirname(getfile(currentframe())))), 'data')

def get_data_dir():
    return _data_dir
