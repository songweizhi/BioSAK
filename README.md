
## BioSAK (A Swiss Army Knife for Biologists)

[![pypi licence       ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version       ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 

Contact
---

+ Weizhi Song
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia
+ E-mail: songwz03@gmail.com

Dependencies
---

+ [BioPython](https://github.com/biopython/biopython.github.io/)

Installation
---

+ BioSAK has been tested on Linux/Mac, but **NOT** supported on Windows.

+ BioSAK is implemented in python3, you can install it with pip3:

      # for the first time installation
      pip3 install BioSAK
      
      # for later updating
      pip3 install --upgrade BioSAK
      
+ For UNSW Katana users

      ############## install BioSAK with Python virtual environment ##############
      
      module load python/3.7.3
      mkdir ~/mypython3env_BioSAK
      python3 -m venv --system-site-packages ~/mypython3env_BioSAK
      source ~/mypython3env_BioSAK/bin/activate
      pip3 install BioSAK
        
      ################################ run BioSAK ################################

      # If you want to run BioSAK later, just run the following commands 
      # to activate the virtual environment and it's ready for running.
      module load python/3.7.3
      source ~/mypython3env_BioSAK/bin/activate
      BioSAK -h

Help information
---

    BioSAK -h
    BioSAK COG2014 -h
    BioSAK select_seq -h
    BioSAK dwnld_GenBank_genome -h

    