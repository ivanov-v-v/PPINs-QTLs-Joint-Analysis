import os
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def ensure_dir(destdir):
    dirptr = os.path.dirname(destdir)
    if not os.path.exists(dirptr):
        os.mkdir(dirptr)

def rgba2hex(rgba_color):
    '''
    Converts rgba to hex.
    Used for making .svg network plots

    :param rgba_color: 4-tuple
    :return: string with hex color representation
    '''
    r = int(rgba_color[0] * 255)
    g = int(rgba_color[1] * 255)
    b = int(rgba_color[2] * 255)
    a = rgba_color[3]

    red = int(((1 - a) * 255) + (a * r))
    green = int(((1 - a) * 255) + (a * g))
    blue = int(((1 - a) * 255) + (a * b))
    return '#{r:02x}{g:02x}{b:02x}'.format(r=red, g=green, b=blue)
