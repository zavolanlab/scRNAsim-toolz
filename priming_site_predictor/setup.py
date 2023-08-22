"""Set up project."""
from pathlib import Path
from setuptools import setup, find_packages

project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as f:
    INSTALL_REQUIRES = f.read().splitlines()

URL = ('https://git.scicore.unibas.ch/zavolan_group/'
       'tools/priming-site-predictor')


setup(
    name='primingsitepredictor',
    url=URL,
    author='Robin Christen & Max Baer',
    author_email='robin.christen@stud.unibas.ch & max.baer@swisstph.ch',
    description='Priming Site Predictor',
    license='MIT,',
    version='0.1.0',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={
        'console_scripts': [
            'primingsitepredictor = primingsitepredictor.cli:main'
        ]
    })
