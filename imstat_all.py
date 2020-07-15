#!/usr/bin/python

import os, string
import numpy as np
import sys

from astropy.io import fits

# this script reads in fits image files and writes output giving statistics 
# about the image

# get fits files
filenames = os.listdir(os.curdir)
for i in range(len(filenames)):
# if '.FITS' in filenames[i] and '-2l_10' in filenames[i]:
  if '.fits' in filenames[i]: 
    fitsimage = filenames[i]
    img = fits.getdata(fitsimage, ext=0)

    mindata=np.amin(img)
    maxdata=np.amax(img)

    minlocation=np.where(img==mindata)
    maxlocation=np.where(img==maxdata)
# Content: minlocation[1] is array of x-coordinates where min occurs
# Content: minlocation[0] is array of y-coordinates where min occurs
# Content: minlocation[1][0] is the x-coordinate of the first occurrence of min
# Content: minlocation[1][1] is the x-coordinate of the second occurrence of min

    sigma=np.std(img)

# If more than one occurrence of min or max, print warning and list the number of occurrences.
    if (len(minlocation[0])>1):
      print('# ',fitsimage,len(minlocation[0]),' pixels at minimum ',mindata)

    if (len(maxlocation[0])>1):
      print('# ',fitsimage,len(maxlocation[0]),' pixels at maximum ',maxdata)

    if ((len(minlocation[0])==1) and (len(maxlocation[0])==1)):
      #print mindata,minlocation[0][0],len(minlocation[0])
      #print maxdata,maxlocation[0][0],
      print("%35s  %10.4e   %15.10e   %5d  %5d    %15.10e   %5d  %5d  %10.4e " % (fitsimage,sigma,mindata,minlocation[1][0],minlocation[0][0],maxdata,maxlocation[1][0],maxlocation[0][0],((maxdata-mindata)/sigma)))


