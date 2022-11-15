from setuptools import setup, find_packages

setup(
    name='sequence_extractor',
    author='Samuel Mondal',
    url='https://git.scicore.unibas.ch/zavolan_group/tools/transcript-sequence-extractor',
    license='MIT',
    version='0.2.0',
    packages=find_packages(),
    install_requires=[],
    entrypoints={
        'console_scripts': ['sequence_extractor=sequence_extractor.cli:main']
        }
)
