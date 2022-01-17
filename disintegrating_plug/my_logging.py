import logging
import os

logger = logging.getLogger('plug')
my_level = os.environ.get('LOGLEVEL', 'INFO').upper()
logger.setLevel(my_level)
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
ch.setLevel(my_level)
logger.addHandler(ch)
