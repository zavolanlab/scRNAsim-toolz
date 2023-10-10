"""scRNAsim-toolz package definition."""
from pathlib import Path
from setuptools import setup, find_packages  # type: ignore
from scRNAsim_toolz.version import __version__

# Read long description from file
with open("README.md", "r", encoding="utf-8") as fh:
    LONG_DESCRIPTION = fh.read()

# Read requirements from file
INSTALL_REQUIRES = []
project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as fh:
    INSTALL_REQUIRES = fh.read().splitlines()

setup(
    name="scrnasim-toolz",
    version=__version__,
    description=(
        "Repository for the tools used by scRNAsim"
    ),
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/zavolanlab/scRNAsim-toolz",
    author="Zavolan Lab, Biozentrum, University of Basel",
    maintainer="Máté Balajti",
    maintainer_email="mate.balajti@unibas.ch",
    project_urls={
        "Repository": "https://github.com/zavolanlab/scRNAsim-toolz",
    },
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    include_package_data=True,
    setup_requires=[
        "setuptools_git == 1.2",
    ],
    entry_points={
        'console_scripts': [
            'transcript-sampler = scRNAsim_toolz.'
            'transcript_sampler.cli:main',
            'structure-generator = scRNAsim_toolz.'
            'structure_generator.cli:app',
            'sequence-extractor = scRNAsim_toolz.'
            'sequence_extractor.cli:main',
            'priming-site-predictor = scRNAsim_toolz.'
            'priming_site_predictor.cli:main',
            'cdna-generator = scRNAsim_toolz.'
            'cdna_generator.cli:main',
            'fragment-selector = scRNAsim_toolz.'
            'fragment_selector.cli:main',
            'read-sequencer = scRNAsim_toolz.'
            'read_sequencer.cli:main',
        ],
    },
)
