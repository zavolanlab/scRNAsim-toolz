"""Set up project."""
from pathlib import Path
from setuptools import setup, find_packages  # type: ignore

project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as f:
    INSTALL_REQUIRES = f.read().splitlines()

URL = ('https://git.scicore.unibas.ch/zavolan_group/'
       'tools/transcript-structure-generator')

setup(
    name='transcript-structure-generator',
    version='0.2.0',
    url=URL,
    license='MIT',
    author='Larissa Glass, Michael Zimmermann, Andri Fraenkl',
    author_email='mate.balajti@unibas.ch',
    description='Transcript structure generator',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={
        'console_scripts': ['transcript-structure-generator=tsg.cli:app']
        }
)
