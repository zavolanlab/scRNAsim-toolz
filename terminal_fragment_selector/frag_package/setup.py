From setup tools import setup, find_packages

setup(
	name = “prog_for_life_science”, 
	url = ‘ssh://git@git.scicore.unibas.ch:2222/zavolan_group/tools/terminal-fragment-selector.git’
	author = ‘Hugo Madge Leon, Sunho Kim, Tanya Nandan’,
	author_email = ‘kimsun@ethz.ch’,
	License = ‘MIT’,
	version = ‘1.0.0’,
	packages = find_packages(),
	install_requires= ['random','argparse 1.4.0','os.path', 'check_positive', 'arghelper 0.4.2', 'fragmentation_v2', 'numpy 1.23.4', 'pandas 1.5']
)