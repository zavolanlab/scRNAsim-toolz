"""Test functions for cli.py."""
import unittest
from unittest.mock import patch
from ...scRNAsim_toolz.priming_site_predictor import cli  # type: ignore
from ...scRNAsim_toolz.priming_site_predictor.cli import (  # type: ignore
    setup_logging
)


class TestCli(unittest.TestCase):
    """Test cli.py methods."""

    def test_parse_args(self):
        """Test parsing --energy-cutoff as a float."""
        # Test parsing a valid float
        float_num = 3.14
        float_str = str(float_num)
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", float_str
                ]):
            args = cli.parse_args()
            self.assertEqual(args.energy_cutoff, float(float_str))

        # Test parsing a float number with leading and trailing whitespace
        float_num = 3.14
        float_str = "  " + str(float_num) + "  "
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", float_str
                ]):
            args = cli.parse_args()
            self.assertEqual(args.energy_cutoff, float(float_str))

        # Test parsing a float number with a decimal point \
        # and no digits after it
        float_num = 3.
        float_str = str(float_num)
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", float_str
                ]):
            args = cli.parse_args()
            self.assertEqual(args.energy_cutoff, float(float_str))

        # Test parsing a float number with a decimal point \
        # and multiple digits after it
        float_num = 3.14159
        float_str = str(float_num)
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", float_str
                ]):
            args = cli.parse_args()
            self.assertEqual(args.energy_cutoff, float(float_str))

        # Test parsing a negative float number
        float_num = -3.14
        float_str = str(float_num)
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", float_str
                ]):
            args = cli.parse_args()
            self.assertEqual(args.energy_cutoff, float(float_str))

    def test_invalid_parse_args(self):
        """Test parsing --energy-cutoff as an invalid float."""
        # Test parsing a float number with a letter in it
        invalid_float_num = "3.14a"
        invalid_float_str = str(invalid_float_num)
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", invalid_float_str
                ]):
            with self.assertRaises(SystemExit):
                cli.parse_args()

        # Test parsing a float number with multiple decimal points
        invalid_float_num = "3.14.159"
        invalid_float_str = str(invalid_float_num)
        with patch("sys.argv", [
            "primingsitepredictor", "--energy-cutoff", invalid_float_str
                ]):
            with self.assertRaises(SystemExit):
                cli.parse_args()


class TestSetupLogging:
    """Test ``setup_logging()`` function."""

    def test_log_level_default(self):
        """Call without args."""
        setup_logging()
