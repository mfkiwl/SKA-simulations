#!/usr/bin/env python

import os
import sys
import numpy
import math 
import csv
import matplotlib.pyplot as plt
import BDA_utility_functions

def main( argv ):

  print(' ')

  Lx, Ly, Lz = BDA_utility_functions.get_ant_pos()     # get antenna x, y, z ITRF positions
  L = numpy.sqrt(Lx*Lx + Ly*Ly + Lz*Lz)
  num_baselines = L.shape[0]

## input parameters for decorrelation calculations
  pi = math.pi
  freq = 700e6                # frequency in MHz
  c = 2.99792e8               # speed of light in m/s
  wavelength = c / freq       # wavelength in m
  dec = numpy.arange(-70,20,20)           # assumed declinations in degrees
  dec = numpy.radians(dec)
  LHAh = numpy.arange(-6,6.1,2)             # local hour angle range in hours
  LHAh = numpy.arange(-6,6.1,0.1)             # local hour angle range in hours
  LHAd = LHAh * 360 / 24;      # local hour angle range in degrees
  LHA = numpy.radians(LHAd)
  
  obs_lon = 21.4439              # meerkat, longitude in degrees
  #obs_lon = 116.659             # askap longitude =  116.659
  
  GMSTd = LHAd - obs_lon         # Greenwich Mean Sidereal Time range in degrees
  GMST = numpy.radians(GMSTd);   
  T = 0.14;                   # correlator dump time in s
  
  l1 = math.sin(math.radians(2.0));                # l-coordinate of the source;
  m1 = math.sin(math.radians(0.0));                # m-coordinate of the source;
  l2 = math.sin(math.radians(0.0));                # l-coordinate of the source;
  m2 = math.sin(math.radians(2.0));                # m-coordinate of the source;
  n1 = math.sqrt(1 - l1*l1 - m1*m1) - 1
  n2 = math.sqrt(1 - l2*l2 - m2*m2) - 1
  print('n1, n2', n1,n2)


  omegaE = 2 * pi / (23 * 3600 + 56 * 60)  # angular velocity of Earth rotation
                                         # in  rad/s

  baseline_zones = [80, 40, 30, 20, 15, 10, 7.5, 5, 3.75, 2.5, 1.875, 1.25, 0.9375, 0.625, 0.5625, 0.375, 0.28125]
  num_zones = len(baseline_zones)
  for i in range(num_zones):
     baseline_zones[i] = baseline_zones[i] * 1000
  baseline_avg = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
  num_avg = baseline_avg[num_zones]
  baseline_table = numpy.zeros((num_baselines,),numpy.int32)
  for i in range(num_baselines):
     baseline_length = L[i]
     num_avg = baseline_avg[num_zones]
     baseline_table[i] = num_avg 
     for j in range(num_zones):
       if baseline_length  > baseline_zones[j]:
         num_avg = baseline_avg[j]
         baseline_table[i] = num_avg 
         break
  baseline_table = baseline_table * baseline_table
## calculate decorrelation assuming natural weighting
  rotationE = omegaE / wavelength
  pi_t = pi * T
  pi_t_const = 100.0 * (pi_t * pi_t /6)
  decor1 = numpy.zeros((len(dec), len(GMST)), numpy.float64);
  decor2 = numpy.zeros((len(dec), len(GMST)), numpy.float64);
  for idx in range(len(dec)):
    print('computing decorrelation for dec ', math.degrees(dec[idx]))
    dec_sin =  math.sin(dec[idx])
    dec_cos =  math.cos(dec[idx])
    for i in range(len(LHA)):
      lha_cos = math.cos(GMST[i])
      lha_sin = math.sin(GMST[i])
      u,v,w,dudt,dvdt,dwdt, u_max, v_max = BDA_utility_functions.uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)

#     decor1[idx,i] = pi_t_const * numpy.mean(baseline_table * numpy.square(dudt * l1 + dvdt * m1));
#     decor2[idx,i] = pi_t_const * numpy.mean(baseline_table * numpy.square(dudt * l2 + dvdt * m2));
      decor1[idx,i] = pi_t_const * numpy.mean(baseline_table * numpy.square(dudt * l1 + dvdt * m1 + dwdt * n1));
      decor2[idx,i] = pi_t_const * numpy.mean(baseline_table * numpy.square(dudt * l2 + dvdt * m2 + dwdt * n2));

    print('computed decorrelation for dec ', math.degrees(dec[idx]))

  
# read in actual data
  print('processing data file ', argv[1])
  dec_values = BDA_utility_functions.getdata(argv[1], argv[2])

  x = LHAh
  plt.xlim(-6, 6)
  colours = ['b', 'r', 'b', 'r', 'g']
  colours = ['g', 'r', 'b', 'r', 'b']
  for i in range(len(dec)):
    plt.plot(x, decor1[i,:], colours[i])
    plt.plot(x, decor2[i,:], colours[i])

# for i in range(len(dec)):
#   plt.plot(x, decor1[i,:])
#   plt.plot(x, decor2[i,:])

  plt.xlabel('Hour Angle (hours)')
  plt.ylabel('Maximum Residual (percent)')
  plt.title('Natural Wt Theoretical Analysis vs SKA Simulator with Scheme 1 BDA')
  ha =  [-4.5, -3, -1.5, 0, 1, 2, 4]
  ecs = ['10', '-10', '-30', '-50', '-70','10', '-10', '-30', '-50', '-70']
  plot_parm = ['bo', 'ro', 'bs', 'r^', 'g^','bo', 'ro', 'bs', 'r^', 'g^']
  for i in range(len(dec_values)):
      x = ha
      y = dec_values[i]
      try:
        if i < 5:
          plt.plot(x, y,plot_parm[i], label= ecs[i])
        else:
          plt.plot(x, y,plot_parm[i])
      except:
        continue
  plt.legend(loc=1)
# plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
  plot_file = 'Figure_4b_plot'
  print('saving figure ', plot_file)
  plt.savefig(plot_file)
  plt.show()
  plt.clf()
  plt.cla()
  plt.close()

#=============================
# argv[1]  incoming ALBUS results file 
if __name__ == "__main__":
  main(sys.argv)
