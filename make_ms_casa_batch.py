#
# Copyright (C) 2008
# ASTRON (Netherlands Foundation for Research in Astronomy)
# and The MeqTree Foundation
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, seg@astron.nl
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# CASA script that uses the CASA simulator to create
# UV tracks for SKA simulations

import random
import numpy
import os
import math
import time
from string import split, strip

sm_a = 6378137.0
invf = 298.257223563
f = 1.0 / invf

# Convert latitude (radians), longitude (radians) and elevation (metres) to ITRF XYZ
def WGS84ToITRF(lat, lon, h): # WGS-84 to ITRF
        SINK = math.sin(lat)
        COSK = math.cos(lat)
        e2 = 2.0 * f - f * f
        v = sm_a / math.sqrt(1.0 - e2 * SINK * SINK)
        x = (v + h) * COSK * math.cos(lon)
        y = (v + h) * COSK * math.sin(lon)
        z = ((1 - e2) * v + h) * SINK
        return x, y, z


def generate_data(msname, DEC0, ha_index):

  print 'generating data for MS ', msname
  # Refer to docs for more details.
  # Feedback on API is welcomed.
  print 'Setting declination to ', DEC0

# basic parameters - things which are most likely to be changed
  observatory  = 'MeerKAT'         # where the telescope is situated
  RA0          = "0h0m1.0"     # field centre RA
# DEC0         = "-50d00m00s"   # field centre DEC
# DEC0         = "-70d00m00s"   # field centre DEC
# DEC0         = "-10d00m00s"   # field centre DEC
  Freq         = "700MHz"     # frequency at lower edge of band
  num_channels = 1             # number of channels
  channel_inc  = '50.0MHz'      # channel increment
  Stokes       = 'XX XY YX YY' # Stokes parameters - other option is 'RR RL LR LL'
  int_time     = '0.14s'         # integration period
  ref_time = me.epoch('UTC','2015/12/12.001')  # a reference time
  starttime = [-16300,-14720,-5500,-10, 3500, 7100,14300] # 1a configuration
  starttime = [-16300,-10900,-5500,-10, 3500, 7100,14300]
  #Seven Hour angles: -4.5h, -3h, -1.5h, 0h, 1h, 2h, 4h.
  #Five Declinations: -70, -50, -30, -10, +10 degrees.

  scanlength   =  72.0      # length of observation, in sec
  scanlength   =  144.0      # length of observation, in sec
  scanlength   =  147.0      # length of observation, in sec
  scanlength   =  420.0      # length of observation, in sec
  noise        =  '0.0Jy'      # add some noise. Use '0.0Jy' if no noise wanted
  FOV          =  5            # field of view in arcmin
  num_sources  = 1             # number of sources to observe, randomized over FOV
  use_meerkat = True           # add meerkat antennas to array
  use_ska1 =  True                # add ska1 antennas to array
# use_ska1 =  False                # add ska1 antennas to array
  
  # first delete old MS, test images, etc
  print '*** deleting previous stuff ***'
  system_call = "/bin/rm -rf " + msname 
  os.system(system_call)
  
  # get telescope positions from text file
  diam=[]
  xx=[]
  yy=[]
  zz=[]
  antnames=[]
  print 'loading antenna positions'
  meerkat_filename = 'meerkat_sa.pos'
  SKA1_filename = 'ska1_mid_sa_BDv2.pos'
  count = 0
  if use_ska1:
     text = open(SKA1_filename, 'r').readlines()
     L = len(text)
     for i in range(L):
#    for i in range(0,L,4):
        try:
           info = split(strip(text[i]))
           if float(info[0]) > 0:
             diam.append(15.0)
             lon = math.radians(float(info[0]))
             lat = math.radians(float(info[1]))
             el = float(info[2])
             x, y, z = WGS84ToITRF(lat, lon, el)
      #      print "  X=%f, Y=%f, Z=%f" %(x, y, z)
      #      print x,y,z
             xx.append(x)
             yy.append(y)
             zz.append(z)
             antnames.append('SKA ' + str(count))
             count = count + 1
        except:
          pass
     print 'using number of SKA antenns = ', len(antnames)
  
  
  if use_meerkat:
     text = open(meerkat_filename, 'r').readlines()
     L = len(text)
     for i in range(L):
#    for i in range(0,L,4):
        try:
           info = split(strip(text[i]))
           if float(info[0]) > 0:
             diam.append(15.0)
             lon = math.radians(float(info[0]))
             lat = math.radians(float(info[1]))
             el = float(info[2])
             x, y, z = WGS84ToITRF(lat, lon, el)
      #      print "  X=%f, Y=%f, Z=%f" %(x, y, z)
      #      print x,y,z
             xx.append(x)
             yy.append(y)
             zz.append(z)
             antnames.append('SKA ' + str(count))
             count = count + 1
        except:
          pass
     print 'imcluding meerkat, number of SKA antenns = ', len(antnames)
  
  sm.open(msname)
  print 'opened MS'
  sm.setspwindow(spwname='SKA', freq=Freq, refcode='LSRK',
		        deltafreq=channel_inc, freqresolution=channel_inc, 
		        nchannels=num_channels, stokes=Stokes)
  pos_obs = me.observatory(observatory)
  sm.setconfig(telescopename=observatory, antname=antnames,
         x=xx, y=yy,z=zz, dishdiameter=diam, mount='ALT-AZ',
         coordsystem='global')
  #      coordsystem='global', referencelocation=pos_obs)
  
  print 'setting up simulator specifications'
  dir0 = me.direction(rf="J2000",  v0=RA0, v1=DEC0)
  sm.setfield(sourcename='SKA_test', sourcedirection=dir0)
  sm.settimes(integrationtime=int_time, usehourangle=True, referencetime=ref_time)
# sm.setlimits(shadowlimit=0.001, elevationlimit='8.0deg')
  sm.setlimits(shadowlimit=0.001, elevationlimit='0.5deg')
  sm.setauto(autocorrwt=0.0)
  
  if ha_index.find('all') > -1:
    start_index = 0
    end_index = len(starttime)
  else:
    start_index = int(ha_index)
    end_index = start_index + 1
  
  for i in range(start_index, end_index):
    print ' **** observing'
    sm.observe('SKA_test', 'SKA', add_observation=True, starttime=str(starttime[i])+'s', 
                                    stoptime=str(starttime[i]+scanlength)+'s')
    me.doframe(ref_time)
    me.doframe(pos_obs)
    sm.setdata(msselect='SCAN_NUMBER=='+ str(i))
  
  if noise != '0.0Jy':
    try:
      sm.setnoise(mode='simplenoise', simplenoise=noise)
      print ' **** corrupting'
      sm.corrupt();
    except:
      print ' **** failure when setting noise for corruption'
      pass
  
  sm.done()
if __name__ == "__main__":
# argv[0] through argv[5] are parameter generated when this script is 
# used with casa_batch
# argv[6] = name of measurement set
# argv[7] = hour angle index
  print '***** ', sys.argv
  try:
    generate_data(sys.argv[6],sys.argv[7],sys.argv[8])
  except:
    print 'trying call for CASA rel 5'
    generate_data(sys.argv[3],sys.argv[4],sys.argv[5])
