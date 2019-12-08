
## BioSAK (A Swiss Army Knife for Biologists)

[![pypi licence       ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version       ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 

Contact
---

+ Weizhi Song (songwz03@gmail.com)
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia

Dependencies
---

+ [BioPython](https://github.com/biopython/biopython.github.io/)

Installation
---

+ BioSAK has been tested on Linux/Mac, but **NOT** supported on Windows.

+ BioSAK is implemented in python3, you can install it with:

      pip3 install BioSAK

+ For UNSW Katana users

      ############## install BioSAK with Python virtual environment ##############
      
      # create a python virtual environment
      module load python/3.7.3
      mkdir ~/mypython3env
      python3 -m venv --system-site-packages ~/mypython3env
      source ~/mypython3env/bin/activate
        
      # for the first time installation
      pip3 install BioSAK

      # for later updating
      pip3 install --upgrade BioSAK

      # to leave Python virtual environment
      deactivate 
        
      ################################ run BioSAK ################################

      # If you want to run BioSAK later, just run the following commands 
      # to activate the virtual environment and BioSAK is ready for running
      module load python/3.7.3
      source ~/mypython3env/bin/activate
      BioSAK -h

Help information
---

    $ BioSAK -h
    $ BioSAK COG2014 -h
    select_seq
