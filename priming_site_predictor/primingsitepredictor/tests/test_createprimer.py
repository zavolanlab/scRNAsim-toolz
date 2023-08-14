
# Imports
import pytest
import createprimer as cp
import unittest
from unittest.mock import mock_open, patch
from examples.count_lines.file_reader import FileReader

# Create an instance of the CreatePrimer class

primer = cp.CreatePrimer( name='primer1', primerlength=15)

# Test for the __init__ method
def test_init(self):
    assert primer.name == primer.name
    assert primer.primer_length == len(primer.primer_length)
    with pytest.raises(ValueError):
        primer.__init__(type(int), type(str))


# Test for the create_fasta method
def test_CreateFasta(self):
    fake_file_path = "fake/file/path"
    content = "Message to write on file to be written"
    with patch('examples.write_on_file.file_writer.open', mock_open()) as mocked_file:
        FileWriter().write(fake_file_path, content)

        # assert if opened file on write mode 'w'
        mocked_file.assert_called_once_with(fake_file_path, 'w')

        # assert if write(content) was called from the file opened
        # in another words, assert if the specific content was written in file
        mocked_file().write.assert_called_once_with(content)