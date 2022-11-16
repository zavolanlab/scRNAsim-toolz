from setuptools import setup
setup(
    name = 'cli-primingsitepredictor',
    url = 'https://git.scicore.unibas.ch/zavolan_group/tools/priming-site-predictor/-/tree/main/CLI',
    author = 'Robin Christen & Max Baer',
    author_email = 'robin.christen@stud.unibas.ch & max.baer@swisstph.ch',
    description = 'Command-Line Interface',
    license = 'MIT,',
    version = '0.1.0',
    packages = ['primingsitepredictor'],
    entry_points = {
        'console_scripts': [
            'primingsitepredictor = primingsitepredictor.cli:letsgo'
        ]
    })
