from setuptools import setup, find_packages
from pathlib import Path
project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt", "r", encoding="utf-8") as _file:
    INSTALL_REQUIRES = _file.read().splitlines()

setup(
    name='readsequencer',
    version='0.1.1',
    url='https://git.scicore.unibas.ch/zavolan_group/tools/read-sequencer',
    license='MIT',
    author='Clara Serger, Michael Sandholzer and Christoph Harmel',
    author_email='christoph.harmel@unibas.ch',
    description='Simulates sequencing with a specified read length from sequences specified by a FASTA file.',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={'console_scripts': ['readsequencer=readsequencer.cli:main']}
)
