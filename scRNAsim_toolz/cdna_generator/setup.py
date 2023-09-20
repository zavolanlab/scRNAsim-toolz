"""Set up project."""
from pathlib import Path
from setuptools import setup, find_packages

project_root_dir = Path(__file__).parent.resolve()

with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as f:
    INSTALL_REQUIRED = f.read().splitlines()

URL = ('https://git.scicore.unibas.ch/zavolan_group/'
       'tools/cdna-generator')

setup(
    name='cdna-generator',
    version='0.1.1',
    url=URL,
    license='MIT',
    author='Eric Boittier, Bastian Wagner, Quentin Badolle',
    author_email='me@email.org',
    description='cDNA generator',
    packages=find_packages(),
    install_required=INSTALL_REQUIRED,
    entry_points={
        'console_scripts': [
            'cdna-generator=cdna.cli:main'
            ]
        }
)
