from setuptools import setup, find_packages

setup(
    name='cdna',
    url='https://gitlab.com/my_user_name/my_package.git',
    author='My Name',
    author_email='me@email.org',
    description='Brief package description',
    license='MIT',
    version='1.0.0',
    packages=find_packages(),  # this will autodetect Python packages from the directory tree, e.g., in `code/`
    install_requires=[],  # add here packages that are required for your package to run, including version or range of versions
)
