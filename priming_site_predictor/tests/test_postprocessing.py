
import unittest
from unittest import mock
import math
from primingsitepredictor.prime_site_predictor import PrimingSitePredictor


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
                'Transcript_1\tRIBlast\tPriming_site\t(0-14:2988-2974)\t+\t.\t+\t.\tInteraction_Energy\t0.75\n'
            )
            print(self.post_processor.generate_gtf())
            print(expected_output)
            # self.assertEqual(self.post_processor.generate_gtf(), expected_output)

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
            "0,Test_Primer,15,Transcript_1,3233,1.49191,-9.76,-8.26809,(0-14:2988-2974)\n",
            "1,Test_Primer,15,Transcript_1,3233,1.02308,-9.76,-8.73692,(0-14:18-4)\n",
        ]

        with mock.patch("builtins.open", mock.mock_open(read_data="")) as m, \
                mock.patch.object(
                    m.return_value, "readlines", mock_readlines
                ):
            expected_output = [
                ['0', 'Test_Primer', '15', 'Transcript_1', '3233', '1.49191', '-9.76', '-8.26809', '(0-14:2988-2974)'],
                ['1', 'Test_Primer', '15', 'Transcript_1', '3233', '1.02308', '-9.76', '-8.73692', '(0-14:18-4)']
            ]
            self.assertEqual(self.post_processor.create_list_from_output(), expected_output)


if __name__ == "__main__":
    unittest.main()

# class TestPrimingSitePredictor(unittest.TestCase):
#     """Test PrimingSitePredictor."""

#     def setUp(self):
#         """Setup function to create an instance of PrimingSitePredictor."""
#         self.post_processor = PrimingSitePredictor()

#     def test_init(self):
#         """Test the __init__() method."""
#         # Test if generate_gtf method is being called
#         with unittest.mock.patch.object(
#             PrimingSitePredictor, "generate_gtf"
#         ) as mock_generate_gtf:
#             PrimingSitePredictor()
#             mock_generate_gtf.assert_called_once()

#         # Test if generate_gtf returns the expected output
#         expected_output = (
#             'Transcript_1\tRIBlast\tPriming_site\t2974[3257 chars]3"\n'
#         )
#         self.assertEqual(self.post_processor.generate_gtf(), expected_output)

#     def test_calculate_energy(self):
#         """Test the calculate_energy() method."""

#         def calculate_energy(value):
#             energy_constant = 1.380649*10**(-23)*298
#             kcalmol_joul = 6.9477*10**-21
#             return (math.exp(-float(value)*kcalmol_joul/energy_constant))

#         # set a decimal place
#         decimalPlace = 20

#         # Test for a positive value
#         self.testinstance = pp.PostProcessRIBlast.calculate_energy(self, 5.0)
#         expected_output = calculate_energy(5.0)
#         self.assertAlmostEqual(self.testinstance, expected_output, decimalPlace)


#         # Test for a negative value
#         self.testinstance = pp.PostProcessRIBlast.calculate_energy(self, -5.0)
#         expected_output = calculate_energy(-5.0)
#         self.assertAlmostEqual(self.testinstance, expected_output, decimalPlace)

#         # Test for zero
#         self.testinstance = pp.PostProcessRIBlast.calculate_energy(self, 0)
#         expected_output = calculate_energy(0)
#         self.assertAlmostEqual(self.testinstance, expected_output, decimalPlace)

#     def test_create_list_from_output(self):
#         """Test the create_list_from_output() method."""
#         # create an instance test of the class PostProcessRIBlast()
#         test = pp.PostProcessRIBlast()
#         test.file = "RIBlast output example.txt"

#         with open('RIBlast output example.txt', 'r', encoding="utf-8") as file:
#             data_list = file.readlines()

#         # Remove the header
#         data_list = data_list[3:]

#         # Convert each row of the list into a list of values
#         expected_output = []
#         for row in data_list:
#             # Split the row by the comma character
#             values = row.strip().split(',')
#             # Append the list of values to the final list
#             expected_output.append(values)
#         return expected_output
#         self.assertEqual(test.create_list_from_output(), expected_output)

#     def test_create_pandas_df(self):
#         """Test the create_pandas_df() method with mock data."""
#         # create a test instance
#         test_instance = pp.PostProcessRIBlast()

#         # create some mock data for the input list
#         test_instance.interaction_list = [
#                                 ['Id', 'Query name', 'Query Length', 'Target name', 'Target Length', 'Accessibility Energy', 'Hybridization Energy', 'Interaction Energy', 'BasePair'],
#                                 [0, 'Test_Primer', 15, 'Transcript_1', 3233, 1.49191, -9.76, -8.26809, '(0-14:2988-2974)'],
#                                 [1, 'Test_Primer', 15, 'Transcript_1', 3233, 1.02308, -9.76, -8.73692, '(0-14:18-4)'],
#                                 [2, 'Test_Primer', 15, 'Transcript_1', 3233, 0.947439, -9.73, -8.78256, '(0-14:17-3)'],
#                                 [3, 'Test_Primer', 15, 'Transcript_1', 3233, 0.793049, -9.73, -8.93695, '(0-14:16-2)'],
#                                 [4, 'Test_Primer', 15, 'Transcript_1', 3233, 0.483869, -9.73, -9.24613, '(0-14:15-1)'],
#                                 [5, 'Test_Primer', 15, 'Transcript_1', 3233, 0.441093, -9.17, -8.72891, '(0-14:14-0)']]

#         # create the DataFrame using the mock data
#         df = test_instance.create_pandas_df()

#         #self.assertEqual(len(df),7)

#         # check that the DataFrame has the expected column names
#         #self.assertEqual(list(df.columns), ['Id','Query name', 'Query Length', 'Target name', 'Target Length', 'Accessibility Energy', 'Hybridization Energy', 'Interaction Energy', 'BasePair'])

#         # check that the values in the 'Id' column are correct
#         self.assertEqual(list(df['Id']), ['Id',0,1,2,3,4,5])

#         # check that the values in the 'Query name' column are correct
#         self.assertEqual(list(df['Query name']), ['Query name','Test_Primer','Test_Primer','Test_Primer','Test_Primer','Test_Primer','Test_Primer'])

#         # check that the values in the 'Query Length' column are correct
#         self.assertEqual(list(df['Query Length']), ['Query Length',15, 15, 15, 15, 15, 15])

#         # check that the values in the 'Target name' column are correct
#         self.assertEqual(list(df['Target name']), ['Target name','Transcript_1','Transcript_1','Transcript_1','Transcript_1','Transcript_1','Transcript_1'])

#         # check that the values in the 'Target Length' column are correct
#         self.assertEqual(list(df['Target Length']), ['Target Length', 3233, 3233, 3233, 3233, 3233, 3233])

#         # check that the values in the 'Accessibility Energy' column are correct
#         self.assertEqual(list(df['Accessibility Energy']), ['Accessibility Energy',1.49191, 1.02308, 0.947439, 0.793049, 0.483869, 0.441093])

#         # check that the values in the 'Hybridization Energy' column are correct
#         self.assertEqual(list(df['Hybridization Energy']), ['Hybridization Energy',-9.76, -9.76, -9.73, -9.73, -9.73, -9.17])

#         # check that the values in the 'Interaction Energy' column are correct
#         self.assertEqual(list(df['Interaction Energy']), ['Interaction Energy',-8.26809, -8.73692, -8.78256, -8.93695, -9.24613, -8.72891])

#         # check that the values in the 'BasePair' column are correct
#         self.assertEqual(list(df['BasePair']), ['BasePair','(0-14:2988-2974)', '(0-14:18-4)', '(0-14:17-3)', '(0-14:16-2)', '(0-14:15-1)', '(0-14:14-0)'])
