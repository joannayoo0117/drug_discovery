from urllib.request import urlretrieve
import os
import logging

DATASET_URL = 'http://deepchem.io.s3-website-us-west-1.amazonaws.com/datasets'
DATASET_TO_FILE = {
    'pdbbind': 'full_smiles_labels.csv',
}
logger = logging.getLogger(__name__)

def download_dataset(save_dir='./data/', type='pdbbind'):
    if type not in DATASET_TO_FILE:
        raise ValueError("Dataset type should be one of: {}".format(
            ", ".join(DATASET_TO_FILE.keys())))

    url = DATASET_URL
    fname = DATASET_TO_FILE[type]

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    if os.path.exists(os.path.join(save_dir, fname)):
        logger.info('Dataset already exists')
    else:
        urlretrieve(os.path.join(url, fname),
                    os.path.join(save_dir, fname))
        logger.info('Successfully downloaded {}'.format(fname))



