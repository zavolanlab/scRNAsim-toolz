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
import postprocessing as pp
import unittest
import math
from unittest import mock


class TestPostProcessRIBlast(unittest.TestCase):
    def test_import(self):
        """Test function for module imports"""
        try:
            import pandas
        except ImportError:
            self.fail("Failed to import pandas")

        try:
            import math
        except ImportError:
            self.fail("Failed to import math")

    def setUp(self):
        """Setup function to create an instance of PostProcessRIBlast()"""

        self.post_processor = pp.PostProcessRIBlast()

    def test_init(self):
        """Test function for the __init__() method"""

        # Test if generate_gtf method is being called
        with unittest.mock.patch.object(
            pp.PostProcessRIBlast, "generate_gtf"
        ) as mock_generate_gtf:
            pp.PostProcessRIBlast()
            mock_generate_gtf.assert_called_once()

        # Test if generate_gtf returns the expected output
        expected_output = (
            'Transcript_1\tRIBlast\tPriming_site\t2974[3257 chars]3"\n'
        )
        self.assertEqual(self.post_processor.generate_gtf(), expected_output)

    def test_calculate_energy(self):
        """Test function for the calculate_energy() method"""

        def calculate_energy(value):
            energy_constant = 1.380649*10**(-23)*298
            kcalmol_joul = 6.9477*10**-21
            return (math.exp(-float(value)*kcalmol_joul/energy_constant))

        # set a decimal place
        decimalPlace = 20

        # Test for a positive value
        self.testinstance = pp.PostProcessRIBlast.calculate_energy(self, 5.0)
        expected_output = calculate_energy(5.0)
        self.assertAlmostEqual(self.testinstance, expected_output, decimalPlace)


        # Test for a negative value
        self.testinstance = pp.PostProcessRIBlast.calculate_energy(self, -5.0)
        expected_output = calculate_energy(-5.0)
        self.assertAlmostEqual(self.testinstance, expected_output, decimalPlace)

        # Test for zero
        self.testinstance = pp.PostProcessRIBlast.calculate_energy(self, 0)
        expected_output = calculate_energy(0)
        self.assertAlmostEqual(self.testinstance, expected_output, decimalPlace)


    def test_create_list_from_output(self):
        """Test function for the create_list_from_output() method"""

        #create an instance test of the class PostProcessRIBlast()
        test = pp.PostProcessRIBlast()
        test.file = "RIBlast output example.txt"

        with open('RIBlast output example.txt', 'r') as file:
            data_list = file.readlines()

        # Remove the header
        data_list = data_list[3:]

        # Convert each row of the list into a list of values
        expected_output = []
        for row in data_list:
            # Split the row by the comma character
            values = row.strip().split(',')
            # Append the list of values to the final list
            expected_output.append(values)
        return expected_output
        self.assertEqual(test.create_list_from_output(), expected_output)


    def test_create_pandas_df(self):
        """Test function for the create_pandas_df() method with mock data"""

        # create a test instance
        test_instance = pp.PostProcessRIBlast()

        # create some mock data for the input list
        test_instance.interaction_list= [
                                ['Id','Query name', 'Query Length', 'Target name', 'Target Length', 'Accessibility Energy', 'Hybridization Energy', 'Interaction Energy', 'BasePair'],
                                [0,'Test_Primer',15,'Transcript_1',3233,1.49191,-9.76,-8.26809,'(0-14:2988-2974)'],
                                [1,'Test_Primer',15,'Transcript_1',3233,1.02308,-9.76,-8.73692,'(0-14:18-4)'],
                                [2,'Test_Primer',15,'Transcript_1',3233,0.947439,-9.73,-8.78256,'(0-14:17-3)'],
                                [3,'Test_Primer',15,'Transcript_1',3233,0.793049,-9.73,-8.93695,'(0-14:16-2)'],
                                [4,'Test_Primer',15,'Transcript_1',3233,0.483869,-9.73,-9.24613,'(0-14:15-1)'],
                                [5,'Test_Primer',15,'Transcript_1',3233,0.441093,-9.17,-8.72891,'(0-14:14-0)']]

        # create the DataFrame using the mock data
        df = test_instance.create_pandas_df()

        #self.assertEqual(len(df),7)

        # check that the DataFrame has the expected column names
        #self.assertEqual(list(df.columns), ['Id','Query name', 'Query Length', 'Target name', 'Target Length', 'Accessibility Energy', 'Hybridization Energy', 'Interaction Energy', 'BasePair'])

        # check that the values in the 'Id' column are correct
        self.assertEqual(list(df['Id']), ['Id',0,1,2,3,4,5])

        # check that the values in the 'Query name' column are correct
        self.assertEqual(list(df['Query name']), ['Query name','Test_Primer','Test_Primer','Test_Primer','Test_Primer','Test_Primer','Test_Primer'])

        # check that the values in the 'Query Length' column are correct
        self.assertEqual(list(df['Query Length']), ['Query Length',15, 15, 15, 15, 15, 15])

        # check that the values in the 'Target name' column are correct
        self.assertEqual(list(df['Target name']), ['Target name','Transcript_1','Transcript_1','Transcript_1','Transcript_1','Transcript_1','Transcript_1'])

        # check that the values in the 'Target Length' column are correct
        self.assertEqual(list(df['Target Length']), ['Target Length', 3233, 3233, 3233, 3233, 3233, 3233])

        # check that the values in the 'Accessibility Energy' column are correct
        self.assertEqual(list(df['Accessibility Energy']), ['Accessibility Energy',1.49191, 1.02308, 0.947439, 0.793049, 0.483869, 0.441093])

        # check that the values in the 'Hybridization Energy' column are correct
        self.assertEqual(list(df['Hybridization Energy']), ['Hybridization Energy',-9.76, -9.76, -9.73, -9.73, -9.73, -9.17])

        # check that the values in the 'Interaction Energy' column are correct
        self.assertEqual(list(df['Interaction Energy']), ['Interaction Energy',-8.26809, -8.73692, -8.78256, -8.93695, -9.24613, -8.72891])

        # check that the values in the 'BasePair' column are correct
        self.assertEqual(list(df['BasePair']), ['BasePair','(0-14:2988-2974)', '(0-14:18-4)', '(0-14:17-3)', '(0-14:16-2)', '(0-14:15-1)', '(0-14:14-0)'])


