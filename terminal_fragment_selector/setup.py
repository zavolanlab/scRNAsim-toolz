"""Set up project."""
from setuptools import setup, find_packages
from pathlib import Path

project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as f:
    INSTALL_REQUIRES = f.read().splitlines()

url = 'https://git.scicore.unibas.ch/zavolan_group/tools/terminal-fragment-selector'

setup(
    name='terminal-fragment-selector',
    version='0.1.1',
    url=url,
    license='MIT',
    author='Hugo Madge Leon, Sunho Kim, Tanya Nandan',
    author_email='hmadge@ethz.ch',
    description='Terminal fragment selector',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES
)
