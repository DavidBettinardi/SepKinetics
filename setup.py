try:
    from setuptools import setup, find_packages
except:
    from disutils.core import setuptools

from codecs import open
from os import path
import sys

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'SepKinetics', '_version.py')) as version_file:
    exec(version_file.read())

with open(path.join(here, 'README.md')) as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

desc = readme + '\n\n' + changelog
try:
    import pypandoc
    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(path.join(here, 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    'numpy',
]

tests_require = [
    'pytest',
    'pytest-cov',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []


setup(
    name='SepKinetics',
    version=__version__,
    description='Deconvolution of UV-vis kinetic spectra',
    author='David J Bettinardi',
    author_email='david.j.bettinardi@gmail.com',
    url='https://github.com/DavidBettinardi/SepKinetics',
    license=MIT,
    python_requires='2.7',
    zip_safe=False,
    packages=find_packages()
    include_package_data= True,

    
)

if __name__ == '__main__':
    setup()
