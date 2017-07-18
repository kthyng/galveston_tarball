'''
Make symbolic links to figures to change order for ffmpeg since backward

NOT NEEDED ANYMORE SINCE I CHANGED TO INDEX-BASED NAMES
'''

import os
from glob import glob

Dirs = glob('figures/*')  # directories

for Dir in Dirs:
    os.chdir(Dir)
    Files = glob('*')  # Files in Dir
    nfiles = len(Files)
    for i, File in enumerate(Files):
        linkname = 'backward-' + str(nfiles-i).zfill(3) + '.png'
        os.symlink(File, linkname)
    os.chdir('../..')
