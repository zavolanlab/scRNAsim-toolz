import readsequencer.cli
import pytest
from cli_test_helpers import ArgvContext, shell
import os
import glob
def test_entrypoint():
    """
    Is entrypoint script installed? (setup.py)
    """
    result = shell('readsequencer --help')
    assert result.exit_code == 0

def test_usage_no_args():
    """
    Does CLI abort w/o arguments, displaying usage instructions?
    """
    with ArgvContext('readsequencer'), pytest.raises(SystemExit):
        readsequencer.cli.main()

    result = shell('readsequencer')

    assert 'usage:' in result.stderr
