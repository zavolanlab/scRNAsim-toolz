from setuptools import setup
setup(
    name = 'cli-primingsitepredictor',
    url = 'https://git.scicore.unibas.ch/zavolan_group/tools/priming-site-predictor/-/tree/main/CLI',
    author = 'Robin_Christen',
    author_email = 'robin.christen@stud.unibas.ch',
    description = 'Command-Line Interface',
    license = 'MIT,',
    version = '0.1.0',
    packages = ['pycli'],
    entry_points = {
        'console_scripts': [
            'pycli = pycli.__main__:main'
        ]
    })
