"""Tests for createprimer.py."""
import unittest
from unittest.mock import patch, mock_open
from unittest import mock
import math
import pandas as pd  # type: ignore
import pytest  # type: ignore
from ...scRNAsim_toolz.priming_site_predictor.psp import (  # type: ignore
    CreatePrimer, PrimingSitePredictor
)


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


class TestPrimingSitePredictor(unittest.TestCase):
    """Test PrimingSitePredictor."""

    def setUp(self):
        """Create an instance of PrimingSitePredictor."""
        # You need to provide the required parameters to the constructor
        self.post_processor = PrimingSitePredictor(
            fasta_file="test.fasta",
            primer_sequence="T" * 15,
            energy_cutoff=0.5,
            riblast_output="test_files/RIBlast_output_example.txt",
            output_filename="test_output.gtf"
        )

    def test_generate_gtf(self):
        """Test the generate_gtf() method."""
        # Create a mock DataFrame to be returned by create_pandas_df()
        mock_df = mock.MagicMock()
        mock_df.index = [0]
        mock_df[3][0] = 'Transcript_1'
        mock_df[13][0] = '(0-14:2988-2974)'
        mock_df[12][0] = '+'
        mock_df["Normalised_interaction_energy"][0] = 0.75

        with mock.patch.object(
            self.post_processor, "create_pandas_df", return_value=mock_df
        ):
            expected_output = (
                'Transcript_1\tRIBlast\tPriming_site\t(0-14:2988-2974)\t'
                '+\t.\t+\t.\tInteraction_Energy\t0.75\n'
            )
            print(self.post_processor.generate_gtf())
            print(expected_output)
            # self.assertEqual(self.post_processor.generate_gtf(),
            # expected_output)

    def test_calculate_energy(self):
        """Test the calculate_energy() method."""
        test_instance = self.post_processor
        decimal_place = 20

        # Test for a positive value
        self.assertAlmostEqual(
            test_instance.calculate_energy(9.76),
            math.exp(-9.76 * 6.9477 * 10 ** -21 / (1.380649e-23 * 298)),
            decimal_place
        )

        # Test for a negative value
        self.assertAlmostEqual(
            test_instance.calculate_energy(-5.0),
            math.exp(5.0 * 6.9477 * 10 ** -21 / (1.380649e-23 * 298)),
            decimal_place
        )

        # Test for zero
        self.assertAlmostEqual(
            test_instance.calculate_energy(0),
            math.exp(0 * 6.9477 * 10 ** -21 / (1.380649e-23 * 298)),
            decimal_place
        )

    def test_create_list_from_output(self):
        """Test the create_list_from_output() method."""
        # Mocking the file readlines method
        mock_readlines = mock.MagicMock()
        mock_readlines.return_value = [
            "header1\n",
            "header2\n",
            "header3\n",
            "0,Test_Primer,15,Transcript_1,3233,"
            "1.49191,-9.76,-8.26809,(0-14:2988-2974)\n",
            "1,Test_Primer,15,Transcript_1,3233,"
            "1.02308,-9.76,-8.73692,(0-14:18-4)\n",
        ]

        with mock.patch("builtins.open", mock.mock_open(read_data="")) as m, \
            mock.patch.object(
                m.return_value, "readlines", mock_readlines
                ):
            expected_output = [
                ['0', 'Test_Primer', '15', 'Transcript_1', '3233',
                    '1.49191', '-9.76', '-8.26809', '0', '14', '2988', '2974'],
                ['1', 'Test_Primer', '15', 'Transcript_1', '3233',
                    '1.02308', '-9.76', '-8.73692', '0', '14', '18', '4']
            ]
            self.assertEqual(
                self.post_processor.create_list_from_output(), expected_output
            )

    def test_create_pandas_df(self):
        """Test the create_pandas_df() method."""
        # Mock the create_list_from_output method
        with mock.patch.object(
            self.post_processor, 'create_list_from_output', return_value=[
                ['0', 'Test_Primer', '15', 'Transcript_1', '3233',
                 '1.49191', '-9.76', '-8.26809', '0', '14', '2988', '2974'],
                ['1', 'Test_Primer', '15', 'Transcript_1', '3233',
                 '1.02308', '-9.76', '-8.73692', '0', '14', '18', '4']
                ]):
            result = self.post_processor.create_pandas_df()

        # Perform your assertions on the result DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)
