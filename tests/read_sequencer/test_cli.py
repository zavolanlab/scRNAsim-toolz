"""Test cli.py."""
import pytest  # type: ignore
from cli_test_helpers import ArgvContext, shell  # type: ignore
from ...scRNAsim_toolz.read_sequencer import cli  # type: ignore


def test_entrypoint():
    """Test if entrypoint script is installed (setup.py)."""
    result = shell('readsequencer --help')
    assert result.exit_code == 0


def test_usage_no_args():
    """Test if CLI aborts w/o arguments, displaying usage instructions."""
    with ArgvContext('readsequencer'), pytest.raises(SystemExit):
        cli.main()

    result = shell('readsequencer')

    assert 'usage:' in result.stderr
