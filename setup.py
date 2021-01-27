import os
from setuptools import setup


def version():

    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'BioSAK', 'VERSION'))
    return version_file.readline().strip()


__long_description__ = ''' BioSAK v%s ''' % version()


setup(name="BioSAK",
      version=version(),
      long_description=__long_description__,
      license="GPL3+",
      author="Weizhi Song",
      author_email="songwz03@gmail.com",
      keywords="Bioinformatics",
      description="BioSAK",
      url="https://github.com/songweizhi/BioSAK",
      packages=['BioSAK'],
      package_data={'': ['*.r', '*.R', '*.py', '*.pl', 'VERSION', '*.hmm']},
      include_package_data=True,
      install_requires=['biopython', 'matplotlib', 'numpy', 'scipy', 'itolapi', 'networkx', 'seaborn', 'lxml', 'beautifulsoup4', 'reportlab', 'ete3'],
      scripts=['bin/BioSAK'])
