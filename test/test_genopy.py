"""This file contains the basic tests for the genopy module
"""
import unittest
import test_settings

# adding the genopy path to PYTHONPATH
##BEGIN appending the path
import sys, os
genopy_path = os.path.dirname(os.getcwd())
sys.path.append(genopy_path)
##END appending the path

class BasicTest(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()