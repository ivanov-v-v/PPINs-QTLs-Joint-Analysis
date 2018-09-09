import os.path

__all__ = ["ROOT_DIR", "LIB_DIR", "BIN_DIR"]

ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
LIB_DIR = os.path.join(ROOT_DIR, "lib")
BIN_DIR = os.path.join(ROOT_DIR, "src")
