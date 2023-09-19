"""Set up project."""
from pathlib import Path
from setuptools import setup, find_packages  # type: ignore

project_root_dir = Path(__file__).parent.resolve()
with open(project_root_dir / "requirements.txt",
          "r", encoding="utf-8") as f:
    INSTALL_REQUIRES = f.read().splitlines()

URL = 'https://git.scicore.unibas.ch/zavolan_group/tools/transcript-sampler'

setup(
    name='transcript-sampler',
    version='0.2.1',
    url=URL,
    license='MIT',
    author='Laura Urbanska, Hugo Gillet, Jakob Rien, Máté Balajti',
    author_email='mate.balajti@unibas.ch',
    description='Transcript sampler',
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    entry_points={
        'console_scripts': ['transcript-sampler=transcript_sampler.cli:main']
        }
)
