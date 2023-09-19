"""Test cli.py."""
import pytest
from cli_test_helpers import ArgvContext, shell  # type: ignore
import readsequencer.cli


def test_entrypoint():
    """Test if entrypoint script is installed (setup.py)."""
    result = shell('readsequencer --help')
    assert result.exit_code == 0


def test_usage_no_args():
    """Test if CLI aborts w/o arguments, displaying usage instructions."""
    with ArgvContext('readsequencer'), pytest.raises(SystemExit):
        readsequencer.cli.main()

    result = shell('readsequencer')

    assert 'usage:' in result.stderr
