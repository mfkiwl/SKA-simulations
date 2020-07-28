#!/usr/bin/env python

import os
import sys
import numpy
import math 
import csv
import matplotlib.pyplot as plt
import BDA_utility_functions

def main( argv ):

  print( ' ')

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


  omegaE = 2 * pi / (23 * 3600 + 56 * 60)  # angular velocity of Earth rotation
                                         # in  rad/s


  baseline_table = numpy.zeros((num_baselines,),numpy.int32)
# Following line give the Stef Salvini algorithm for determining the number of points to average
  eps = numpy.spacing(1)
  maxbmult = 500  
  maxbmult = 48
  maxbmult = 250  
  maxbmult = 125  
  maxbmult = 63  
  maxbmult = 32  
  maxbs = [500, 250, 125, 62, 31]
  maxbs = [500, 250]
  maxbs = [250]
  maxbs = [32]
  for j in range(len(maxbs)):
    maxbmult = maxbs[j]
    print( 'using maxbmult = ', maxbmult)
    bmax = L.max()
    print( 'max baseline ', bmax)
    for i in range(num_baselines):
      baseline_length = L[i]
      num_avg = min(maxbmult,int(numpy.fix((bmax+100*eps)/baseline_length)))
      baseline_table[i] = num_avg

    baseline_table = baseline_table * baseline_table
## calculate decorrelation assuming natural weighting
    rotationE = omegaE / wavelength
    pi_t = pi * T
    pi_t_const = 100.0 * (pi_t * pi_t /6)
    decor1 = numpy.zeros((len(dec), len(GMST)), numpy.float64);
    decor2 = numpy.zeros((len(dec), len(GMST)), numpy.float64);
    for idx in range(len(dec)):
      dec_sin =  math.sin(dec[idx])
      dec_cos =  math.cos(dec[idx])
      for j in range(len(GMST)):
        lha_cos = math.cos(GMST[j])
        lha_sin = math.sin(GMST[j])
        u,v,w,dudt,dvdt,dwdt, u_max, v_max = BDA_utility_functions.uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)

        decor1[idx,j] = pi_t_const * numpy.mean(baseline_table * numpy.square(dudt * l1 + dvdt * m1 + dwdt * n1));
        decor2[idx,j] = pi_t_const * numpy.mean(baseline_table * numpy.square(dudt * l2 + dvdt * m2 + dwdt * n2));

    print( 'processing data file ', argv[1])
    dec_values = BDA_utility_functions.getdata(argv[1], argv[2])

    x = LHAh
    plt.xlim(-6, 6)
#   plt.ylim(0, 1.2)
    colours = ['b', 'r', 'b', 'r', 'g']
    colours = ['g', 'r', 'b', 'r', 'b']
    for i in range(len(dec)):
      plt.plot(x, decor1[i,:], colours[i])
      plt.plot(x, decor2[i,:], colours[i])

#   for i in range(len(dec)):
#     plt.plot(x, decor1[i,:])
#     plt.plot(x, decor2[i,:])

    plt.xlabel('Hour Angle (hours)')
    plt.ylabel('(Image Range) / sigma')
    plt.ylabel('Maximum Difference (percent)')
    plt.title('Natural Wt Theoretical Analysis vs SKA Simulator with Scheme 2 BDA'+ str(maxbmult))
    
    ha =  [-4.5, -3, -1.5, 0, 1, 2, 4]
    decs = ['10','-10', '-30', '-50','-70']
    decs = ['-70', '-50', '-30', '-10', '10']
    ecs = ['10', '-10', '-30', '-50', '-70','10', '-10', '-30', '-50', '-70']
    plot_parm = ['ro', 'g--', 'bs', 'r^', 'g']
    plot_parm = ['bo', 'ro', 'bs', 'r^', 'g']
    plot_parm = ['bo', 'ro', 'bs', 'r^', 'g^','bo', 'ro', 'bs', 'r^', 'g^']
    for i in range(len(dec_values)):
      x = ha
      y = dec_values[i] 
      if i < 5:
        plt.plot(x, y,plot_parm[i], label= ecs[i])
      else:
        plt.plot(x, y,plot_parm[i])
    plot_file = 'Figure_9a_plot'
#   plt.legend(loc=1)
    plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
    print( 'saving figure ', plot_file)
    plt.savefig(plot_file)
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()
#=============================
# argv[1]  incoming ALBUS results file 
if __name__ == "__main__":
  main(sys.argv)
