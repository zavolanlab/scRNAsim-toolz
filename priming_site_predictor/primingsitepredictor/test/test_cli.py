"""
Created on Tue Dec 20 14:06:20 2022

@author: RobinC
"""

# Imports
import os, sys
import pytest
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import createprimer as cp
import unittest
import cli
import argparse
import pytest
from unittest.mock import patch
from unittest.mock import Mock
import io
import sys

class TestCli(unittest.TestCase):

    def test_import(self):
        """Test function for module imports"""
        try:
            import argparse
        except ImportError:
            self.fail("Failed to import argparse")

        try:
            import logging
        except ImportError:
            self.fail("Failed to import logging")

        try:
            import main
        except ImportError:
            self.fail("Failed to import main")



    def test_create_parser(self):

        # Test parsing a float number
        float_num = 3.14
        float_str = str(float_num)
        self.assertEqual(float_num, cli.CLI.create_parser(['--float', float_str]))

        # Test parsing a float number with leading and trailing whitespace
        float_num = 3.14
        float_str = "  " + str(float_num) + "  "
        self.assertEqual(float_num, cli.CLI.create_parser(['--float', float_str]))

        # Test parsing a float number with a decimal point and no digits after it
        float_num = 3.
        float_str = str(float_num)
        self.assertEqual(float_num, cli.CLI.create_parser(['--float', float_str]))

        # Test parsing a float number with a decimal point and multiple digits after it
        float_num = 3.14159
        float_str = str(float_num)
        self.assertEqual(float_num, cli.CLI.create_parser(['--float', float_str]))

        # Test parsing a negative float number
        float_num = -3.14
        float_str = str(float_num)
        self.assertEqual(float_num, cli.CLI.create_parser(['--float', float_str]))

    def test_invalid_create_parser(self):

        # Test parsing a float number with a letter in it
        args = ['--float', '3.14a']
        with self.assertRaises(SystemExit):
            cli.CLI.create_parser(args)

        # Test parsing a float number with multiple decimal points
        args = ['--float', '3.14.159']
        with self.assertRaises(SystemExit):
            cli.CLI.create_parser(args)

    def test_letsgo(self):

        # Test printing a float number
        float_num = 3.14
        float_str = str(float_num)
        with self.assertLogs('', level='INFO') as cm:
            cli.CLI.letsgo(['--float', float_str])
        self.assertEqual(cm.output, [f'INFO:root:Your energy cutoff is {float_num}'])













