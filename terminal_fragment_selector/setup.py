"""Set up project."""
from pathlib import Path
from setuptools import setup, find_packages

project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as f:
    INSTALL_REQUIRES = f.read().splitlines()

URL = ('https://git.scicore.unibas.ch/zavolan_group/'
       'tools/terminal-fragment-selector')

setup(
    name='terminal-fragment-selector',
    version='0.1.1',
    url=URL,
    license='MIT',
    author='Hugo Madge Leon, Sunho Kim, Tanya Nandan',
    author_email='hmadge@ethz.ch',
    description='Terminal fragment selector',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={
        'console_scripts': [
            'terminal-fragment-selector=term_frag_sel.cli:main'
            ]
        }
)
