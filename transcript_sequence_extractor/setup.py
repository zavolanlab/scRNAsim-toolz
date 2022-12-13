from setuptools import setup, find_packages
from pathlib import Path

project_root_dir = Path(__file__).parent.resolve()

with open(project_root_dir / "requirements.txt", "r", encoding="utf-8") as _file:
    INSTALL_REQUIRES = _file.read().splitlines()

setup(
    name='sequence_extractor',
    author='Samuel Mondal',
    author_email='samuel.mondal@unibas.ch',
    url='https://git.scicore.unibas.ch/zavolan_group/tools/transcript-sequence-extractor',
    license='MIT',
    version='0.0.1',
    description='Extracts transcript sequences from gtf file and adds polyA tail to the output sequence',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entrypoints={
        'console_scripts': ['sequence_extractor=sequence_extractor.cli:main']
        }
)
