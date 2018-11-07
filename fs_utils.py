import os.path

__all__ = ["ROOT_DIR", "LIB_DIR", "SRC_DIR", "DATA_DIR"]

ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

LIB_DIR = os.path.join(ROOT_DIR, "lib")
SRC_DIR = os.path.join(ROOT_DIR, "src")
DATA_DIR = os.path.join(ROOT_DIR, "data")
