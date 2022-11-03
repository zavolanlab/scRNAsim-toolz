"""Setup the package."""

from setuptools import setup

setup(
    name='tsg',
    author='Zimmermann, M; Fraenkl, A;Glass, L',
    url='https://git.scicore.unibas.ch/zavolan_group/tools/transcript-structure-generator',
    license='MIT',
    version='0.0.1',
    packages=['tsg'],
    install_requires=['pandas'],
    entry_points={
        'console_scripts': [
            'transcript-generator = tsg:cli',
        ]
    }
)
