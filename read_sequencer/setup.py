 from setuptools import setup

 setup(
   name='awesome_read_sequencer',
   version='0.1.0',
   author='An Awesome Coder',
   author_email='aac@example.com',
   packages=['random', 'sys'],
   scripts=['read_in_FASTA.py','read_sequence.py']
   license='LICENSE.txt',
   description='An awesome package that simulates sequencing of a FASTA file.',
   long_description=open('README.md').read(),
   install_requires=[
       "random",
       "sys"
   ],
)