"""Set up project."""
from pathlib import Path
from setuptools import setup, find_packages  # type: ignore

project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as file:
    INSTALL_REQUIRES = file.read().splitlines()

URL = ('https://git.scicore.unibas.ch/zavolan_group/'
       'tools/transcript-sequence-extractor')

setup(
    name='transcript-sequence-extractor',
    version='0.1.0',
    url=URL,
    license='MIT',
    author='Samuel Mondal',
    author_email='samuel.mondal@unibas.ch',
    description=('Extracts transcript sequences from gtf file'
                 'and adds polyA tail to the output sequence'),
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={
        'console_scripts': [
            'sequence-extractor=sequence_extractor.cli:main'
            ]
        }
)
