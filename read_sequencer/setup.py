from setuptools import setup

setup(
    name='awesome_read_sequencer',
    version='0.1.0',
    author='An Awesome Coder',
    author_email='aac@example.com',
    packages=['random'],
    scripts=['cli.py', 'modules.py'],
    license='LICENSE.txt',
    description='An awesome package that simulates sequencing from sequences specified by a FASTA file.',
    long_description=open('README.md').read(),
    install_requires=['random', 'sys'],
    entry_points={
        'console_scripts': ['read_sequencer=read_sequencer_package/cli.py:main']
    }
)
