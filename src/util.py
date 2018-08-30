import os
import sys

import imageio
import numpy as np
import scipy as sp

roman2arabic = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
    "V": 5,
    "VI": 6,
    "VII": 7,
    "VIII": 8,
    "IX": 9,
    "X": 10,
    "XI": 11,
    "XII": 12,
    "XIII": 13,
    "XIV": 14,
    "XV": 15,
    "XVI": 16
}
arabic2roman = {val: key for key, val in roman2arabic.items()}


def as_roman(chromosome_id):
    """ Warning: yeast-specific function """
    return arabic2roman[chromosome_id]


def as_arabic(chromosome_id):
    """ Warning: yeast-specific function """
    return roman2arabic[chromosome_id]


def eprint(*args, **kwargs):
    """ Print to stderr """
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(destdir):
    """ Create a directory if it doesn't exist """
    if not os.path.exists(destdir):
        os.makedirs(destdir)


def try_imread(filepath, mode='RGBA'):
    """ Read an image if it exists """
    try:
        img = imageio.imread(uri=filepath, format="png", pilmode=mode)
    except FileNotFoundError:
        img = None
    return img


def rgba2hex(rgba_color):
    """
    Converts rgba to hex.
    Used for making .svg network plots

    :param rgba_color: 4-tuple
    :return: string with hex color representation
    """
    r = int(rgba_color[0] * 255)
    g = int(rgba_color[1] * 255)
    b = int(rgba_color[2] * 255)
    a = rgba_color[3]

    red = int(((1 - a) * 255) + (a * r))
    green = int(((1 - a) * 255) + (a * g))
    blue = int(((1 - a) * 255) + (a * b))
    return '#{r:02x}{g:02x}{b:02x}'.format(r=red, g=green, b=blue)


def confidence_interval(data, quantile=0.95):
    """ https://kite.com/python/examples/702/scipy-compute-a-confidence-interval-from-a-dataset """
    n = len(data)
    m = np.mean(data)
    std_err = sp.stats.sem(data)
    h = std_err * sp.stats.t.ppf((1 + quantile) / 2, n - 1)
    return (m - h, m + h)


def bincount_scott(data):
    """ Calculate optimal number of histogram bins given data """
    n = len(data)
    std_err = np.std(data)
    data_range = np.max(data) - np.min(data)
    return int(np.ceil(data_range * (n ** (1 / 3)) / (3.49 * std_err)))


def sample_combinations(dims, nsamp):
    """ TODO """
    idx = np.random.RandomState().choice(np.prod(dims), nsamp, replace=False)
    return np.vstack(np.unravel_index(idx, dims)).T