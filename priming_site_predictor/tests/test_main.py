"""
Created on Tue Dec 20 14:06:20 2022

@author: RobinC
"""

# Imports
import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import unittest
import main as mn
import createprimer as cp
import postprocessing as pp
from io import StringIO

class TestMain(unittest.TestCase):

    def test_import(self):
        """Test function for module imports"""
        try:
            import CreatePrimer
        except ImportError:
            self.fail("Failed to import CreatePrimer")

        try:
            import PostProcessRIBlast
        except ImportError:
            self.fail("Failed to import PostProcessRIBlast")

    def test_generate_riblast_input(self):
        """Test funciton for the method generate_riblast_input()"""

        # Call the CreatePrimer() method and the create_fasta() method
        primer = cp.CreatePrimer()
        primer.create_fasta()

        # get the name and the transcript filename
        primer_filename = primer.name+".fasta"
        transcript_filename = "transcripts.fasta"

        # Call the generate_riblast_input() method
        result = mn.generate_riblast_input()

        assert result == [primer_filename, transcript_filename]


    def test_create_gtf(self):
        """Test funciton for the method create_gtf()"""

        # Call the output of the PostprocessRIBlast() method
        result = pp.PostProcessRIBlast().output

        # Open the expected output
        with open('output_transcripts_df.txt', 'r') as f:
            expected_output = f.read()

        # Use an assertion to check the output result
        assert result == expected_output

        # Capture the output of the print statement
        captured_output = StringIO()
        sys.stdout = captured_output

        # Call the function that contains the print statement
        mn.create_gtf()

        # Check the captured output against the expected value
        assert captured_output.getvalue().strip() == expected_output


    def test_main(self):
        """Test function for the method main()"""

        # Call the CreatePrimer() method and the create_fasta() method
        primer = cp.CreatePrimer()
        primer.create_fasta()

        # get the name and the transcript filename
        primer_filename = primer.name+".fasta"
        transcripts_filename = "transcripts.fasta"

        # Get the results for the called methods
        result1 = mn.generate_riblast_input()
        result2 = mn.create_gtf()

        # Open the expected output
        with open('output_transcripts_df.txt', 'r') as f:
            expected_output = f.read()

        assert result1 == [primer_filename, transcripts_filename]
        assert result2 == expected_output
