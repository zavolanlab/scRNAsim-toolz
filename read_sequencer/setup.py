from setuptools import setup, find_packages

setup(
    name='read_sequencer',
    version='0.1.1',
    url='https://git.scicore.unibas.ch/zavolan_group/tools/read-sequencer',
    license='MIT',
    author='Clara Serger, Michael Sandholzer and Christoph Harmel',
    author_email='christoph.harmel@unibas.ch',
    description='Simulates sequencing with a specified read length from sequences specified by a FASTA file.',
    packages=find_packages(),
    install_requires=['random','textwrap','argparse','logging'],
    entry_points={'console_scripts': ['read_sequencer=read_sequencer_package.cli:main']}
)
