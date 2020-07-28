#!/usr/bin/env python

import os
import time
import sys
import numpy
import math 
import csv
import matplotlib.pyplot as plt
import BDA_utility_functions

def main( argv ):

# set the following variable to True for running under one processor
  single_proc = False
  
  print( ' ')
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  process_start = time.time()
  print( "Start at %s" % startime)

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

  baseline_table = numpy.zeros((num_baselines,),numpy.int32)

# Following line give the Stef Salvini algorithm for determining the number of points to average
  eps = numpy.spacing(1)
  maxbs = [500, 250, 125, 62, 32]
  maxbs = [32]
  maxbs = [250]
  for j in range(len(maxbs)):
    maxbmult = maxbs[j]
    print( 'using maxbmult = ', maxbmult)
    bmax = L.max() / 1000.0
    print( 'max baseline ', bmax)
    for i in range(num_baselines):
      baseline_length = L[i] / 1000.0
      num_avg = min(maxbmult,int(numpy.fix((bmax+100*eps)/baseline_length)))
      baseline_table[i] = num_avg

    baseline_num = baseline_table
    baseline_table = baseline_table * baseline_table
    print( 'baseline_num ', baseline_num)

    threshhold = -1
    threshhold = 20
    w_factor_weight = 0.15
    if single_proc:
      weight_uni = BDA_utility_functions.calculate_basic_uniform_BDA_weight_single(Lx,Ly,Lz,dec,GMST,wavelength, threshhold)
    else:
       weight_uni, weight_uni1 = BDA_utility_functions.calculate_basic_uniform_BDA_weight(Lx,Ly,Lz,dec,GMST,wavelength,threshhold,w_factor_weight)

## calculate basic decorrelation assuming uniform weighting
    rotationE = omegaE / wavelength
    pi_t = pi * T
    pi_t_const = 100.0 * (pi_t * pi_t /6)

# calculate decorrelation with BDA assuming uniform weighting
    decorBDA1uni = numpy.zeros((len(dec), len(GMST)), numpy.float64);
    decorBDA2uni = numpy.zeros((len(dec), len(GMST)), numpy.float64);
    for idx in range(len(dec)):
      dec_sin =  math.sin(dec[idx])
      dec_cos =  math.cos(dec[idx])
      for j in range(len(GMST)):
        lha_cos = math.cos(GMST[j])
        lha_sin = math.sin(GMST[j])
        u,v,w,dudt,dvdt,dwdt, u_max, v_max = BDA_utility_functions.uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)

        weight = numpy.squeeze(weight_uni[idx,j,:])
        weight = weight / baseline_num 
        sum_wt = numpy.sum(weight)
#       inv_sum = pi_t_const / sum_wt 
#       # calculate decorrelation
#       temp1 = baseline_table * numpy.square(dudt * l1 + dvdt * m1 + dwdt * n1);
#       temp2 = baseline_table * numpy.square(dudt * l2 + dvdt * m2 + dwdt * n2);
#       decorBDA1uni[idx, j] = inv_sum * numpy.sum(weight * temp1)
#       decorBDA2uni[idx, j] = inv_sum * numpy.sum(weight * temp2)

# sin(x)/x test
        temp1 = pi_t * baseline_num * (dudt * l1 + dvdt * m1 + dwdt * n1)
        temp2 = pi_t * baseline_num * (dudt * l2 + dvdt * m2 + dwdt * n2);
        temp1 = numpy.sin(temp1) /temp1
        temp2 = numpy.sin(temp2) /temp2
#       decorBDA1uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp1)/sum_wt) - decorBDA1uni[idx, j]  # differences are at most 0.003 percent!
#       decorBDA2uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp2)/sum_wt) - decorBDA2uni[idx, j]  # differences are at most 0.003 percent!
        decorBDA1uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp1)/sum_wt)
        decorBDA2uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp2)/sum_wt)


    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print( "End at %s" % endtime)
    process_end = time.time()
    duration = (process_end - process_start)/3600.0
    print( "pre-plot Total run time: %7.2f hours" % duration)

#   print( 'processing data file ', argv[1])
    dec_values = BDA_utility_functions.getdata(argv[1], argv[2])
#   print( 'observed dec_values ', dec_values)
    ha =  [-4.5, -3, -1.5, 0, 1, 2, 4]
    ecs = ['10', '-10', '-30', '-50', '-70','10', '-10', '-30', '-50', '-70']
    plt.xlim(-6, 6)
    plot_parm = ['ro', 'g--', 'bs', 'r^', 'g']
    plot_parm = ['bo', 'ro', 'bs', 'r^', 'g']
    plot_parm = ['bo', 'ro', 'bs', 'r^', 'g^','bo', 'ro', 'bs', 'r^', 'g^']
    plt.title('Uniform Wt Theoretical Analysis vs SKA Simulator with Scheme 2 BDA'+ str(maxbmult))
    plt.xlabel('Hour Angle (hours) ')
    plt.ylabel('Maximum Difference (percent)')
    x = LHAh
    colours = ['b', 'r', 'b', 'r', 'g']
    colours = ['g', 'r', 'b', 'r', 'b']
    for i in range(len(dec)):
      plt.plot(x, decorBDA1uni[i,:], colours[i])
      plt.plot(x, decorBDA2uni[i,:], colours[i])

#   for i in range(len(dec)):
#     plt.plot(x, decorBDA1uni[i,:])
#     plt.plot(x, decorBDA2uni[i,:])

    for i in range(len(dec_values)):
      x = ha
      y = dec_values[i] 
      if i < 5:
        plt.plot(x, y,plot_parm[i], label= ecs[i])
      else:
        plt.plot(x, y,plot_parm[i])
#   plt.legend(loc=1)
    plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
    plot_file = 'salvini_uniform_weighting_BDA' + str(maxbmult) + '_-1'
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
