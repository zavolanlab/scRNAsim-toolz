"""Tests for createprimer.py."""
import unittest
from unittest.mock import patch, mock_open
import pytest
from primingsitepredictor.prime_site_predictor import CreatePrimer


class TestCreatePrimer(unittest.TestCase):
    """Test CreatePrimer function."""

    def test_init(self):
        """Test the  __init__() default values."""
        primer = CreatePrimer()
        self.assertEqual(primer.name, 'primer')
        self.assertEqual(primer.primer_length, 15)
        self.assertEqual(primer.primer_sequence, 'T'*15)
        self.assertEqual(primer.lines, ['<primer', 'TTTTTTTTTTTTTTT'])
        with pytest.raises(TypeError):
            CreatePrimer(type(int), type(str))

    def test_init_custom(self):
        """Test the  __init__() custom values."""
        primer = CreatePrimer('my_primer', 20)
        self.assertEqual(primer.name, 'my_primer')
        self.assertEqual(primer.primer_length, 20)
        self.assertEqual(primer.primer_sequence, 'T'*20)
        self.assertEqual(primer.lines, ['<my_primer', 'TTTTTTTTTTTTTTTTTTTT'])

    def test_create_fasta(self):
        """Test the test_create_fasta default values."""
        with patch('builtins.open', new_callable=mock_open) as mock_file:
            primer = CreatePrimer()
            primer.create_fasta()

            # assert if opened file on write mode 'w'
            mock_file.assert_called_once_with(
                'primer.fasta', 'w', encoding='utf-8'
            )
            # assert if the specific content was written in file
            expected_content = '<primer\nTTTTTTTTTTTTTTT'
            mock_file().write.assert_called_once_with(expected_content)

    def test_create_fasta_custom(self):
        """Test the test_create_fasta custom values."""
        with patch('builtins.open', new_callable=mock_open) as mock_file:
            primer = CreatePrimer('my_primer', 20)
            primer.create_fasta()

            # assert if opened file on write mode 'w'
            mock_file.assert_called_once_with(
                'my_primer.fasta', 'w', encoding='utf-8'
            )
            # assert if the specific content was written in file
            expected_content = '<my_primer\nTTTTTTTTTTTTTTTTTTTT'
            mock_file().write.assert_called_once_with(expected_content)
