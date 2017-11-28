#!/usr/bin/env python

import os
import urllib

out = 'data/'
files = ['galaxies_DR9_CMASS_North.fits',
         'galaxies_DR9_CMASS_South.fits',
         'randoms_DR9_CMASS_North.fits',
         'randoms_DR9_CMASS_South.fits']
url = "https://data.sdss.org/sas/dr9/boss/lss/"

if not os.path.exists(out):
    os.makedirs(out)
    
for f in files:
    print 'Downloading', url+f
    urllib.urlretrieve (url+f, out+f)
    print "to", out+f
    print
