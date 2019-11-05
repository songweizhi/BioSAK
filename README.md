
## BioSAK 

[![pypi   licence        ](https://img.shields.io/pypi/l/BioSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi   version        ](https://img.shields.io/pypi/v/BioSAK.svg)](https://pypi.python.org/pypi/BioSAK) 


Contact
---

+ Weizhi Song (songwz03@gmail.com)
+ The Centre for Marine Bio-Innovation, University of New South Wales, Sydney, Australia


Dependencies
---

+ [BioPython](https://github.com/biopython/biopython.github.io/)


Installation
---

+ BioSAK is implemented in python3, you can install it with:

        pip3 install BioSAK


+ For UNSW Katana users

        ############################## install BioSAK ##############################
        
        # create a python virtual environment
        module load python/3.7.3
        mkdir ~/mypython3env
        python3 -m venv --system-site-packages ~/mypython3env
        source ~/mypython3env/bin/activate
        
        # to install 
        pip3 install BioSAK
        
        # to upgrade
        pip3 install --upgrade BioSAK
        
        ################################ run BioSAK ################################

        # activate python virtual environment
        module load python/3.7.3
        source ~/mypython3env/bin/activate
        
        # BioSAK is ready to run now
        BioSAK -h


Help information
---

    $ BioSAK -h
