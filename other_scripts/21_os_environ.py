import os
from pathlib import Path


db_var = "ENRICHM_DB"
if db_var in os.environ:
    DATABASE_DIR = os.environ[db_var]
else:
    DATA_PATH = str(Path.home())
    DATABASE_DIR = os.path.join(DATA_PATH, 'databases')

