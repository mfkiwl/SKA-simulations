## Prediction for the decorrelation caused by time averaging
#
# This script predicts the decorrelation caused by time averaging for the
# SKA-mid configuration. This prediction is made for the case with and
# without baseline dependent averaging (BDA).
#
# SJW, 20 January 2017

import random
import numpy
import os
import sys
import math
import time
import matplotlib.pyplot as plt
import BDA_utility_functions


def main(argv):

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
  dec = [-70]
  dec = numpy.radians(dec)
  LHAh = numpy.arange(-6,6.1,2)             # local hour angle range in hours with 2 hr increment
  LHAh = numpy.arange(-6,6.1,0.1)           # local hour angle range in hours with 0.1 hr increment
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
  n1 = math.sqrt(1 - l1*l1 - m1*m1) -1
  n2 = math.sqrt(1 - l2*l2 - m2*m2) -1

  omegaE = 2 * pi / (23 * 3600 + 56 * 60)  # angular velocity of Earth rotation
                                         # in  rad/s

  baseline_zones = [80, 40, 30, 20, 15, 10, 7.5, 5, 3.75, 2.5, 1.875, 1.25, 0.9375, 0.625, 0.5625, 0.375, 0.28125]
  num_zones = len(baseline_zones)
  for i in range(num_zones):
     baseline_zones[i] = baseline_zones[i] * 1000
  baseline_avg = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
  num_avg = baseline_avg[num_zones]
  baseline_table = numpy.zeros((num_baselines,),numpy.int32)
  num_bins = len(baseline_avg)
  num_avg = baseline_avg[num_zones]
  print( 'max avg ', num_avg)
  
  for i in range(num_baselines):
     baseline_length = L[i]
     num_avg = baseline_avg[num_zones]
     baseline_table[i] = num_avg 
     for j in range(num_zones):
       if baseline_length  > baseline_zones[j]:
         num_avg = baseline_avg[j]
         baseline_table[i] =  num_avg 
         break
  baseline_num = baseline_table
  baseline_table = baseline_table * baseline_table

  if single_proc:
    weight_uni = BDA_utility_functions.calculate_basic_uniform_BDA_weight_single(Lx,Ly,Lz,dec,GMST,wavelength)
  else:
     weight_uni, weight_uni1 = BDA_utility_functions.calculate_basic_uniform_BDA_weight(Lx,Ly,Lz,dec,GMST,wavelength)

## calculate basic decorrelation assuming uniform weighting
  rotationE = omegaE / wavelength
  pi_t = pi * T
  pi_t_const = 100.0 * (pi_t * pi_t /6)

# sinc_factor = baseline_num * pi_t  / freq
# print( 'sinc_factor ', sinc_factor
# atten_factor = numpy.sin(sinc_factor) / sinc_factor
# print( 'atten_factor ', atten_factor

# decorrelation with BDA assuming uniform weighting for Willis algorithm
  decorBDA1uni = numpy.zeros((len(dec), len(GMST)), numpy.float64);
  decorBDA2uni = numpy.zeros((len(dec), len(GMST)), numpy.float64);
  for idx in range(len(dec)):
    dec_sin =  math.sin(dec[idx])
    dec_cos =  math.cos(dec[idx])
    for j in range(len(GMST)):
      lha_cos = math.cos(GMST[j])
      lha_sin = math.sin(GMST[j])
      # weights approximuated based on radial density of (u, v)-points
      u,v,w,dudt,dvdt,dwdt,u_max,v_max = BDA_utility_functions.uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)
      weight = weight_uni[idx,j,:]
      weight = weight / baseline_num 
      sum_wt = numpy.sum(weight)
      inv_sum = pi_t_const / sum_wt
    # calculate decorrelation
#     temp1 = baseline_table * numpy.square(dudt * l1 + dvdt * m1);
#     temp2 = baseline_table * numpy.square(dudt * l2 + dvdt * m2);
#     temp1 = baseline_table * numpy.square(dudt * l1 + dvdt * m1 + dwdt * n1);
#     temp2 = baseline_table * numpy.square(dudt * l2 + dvdt * m2 + dwdt * n2);
#     decorBDA1uni[idx, j] = inv_sum * numpy.sum(weight * temp1) 
#     decorBDA2uni[idx, j] = inv_sum * numpy.sum(weight * temp2)

# sin(x)/x test
      temp1 = pi_t * baseline_num * (dudt * l1 + dvdt * m1 + dwdt * n1)
      temp2 = pi_t * baseline_num * (dudt * l2 + dvdt * m2 + dwdt * n2);
#     temp1 = pi_t * baseline_num * (dudt * l1 + dvdt * m1 ) # dwdt does have some impact
#     temp2 = pi_t * baseline_num * (dudt * l2 + dvdt * m2 );
      temp1 = numpy.sin(temp1) /temp1
      temp2 = numpy.sin(temp2) /temp2
#     decorBDA1uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp1)/sum_wt) - decorBDA1uni[idx, j]  # differences are at most 0.003 percent!
#     decorBDA2uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp2)/sum_wt) - decorBDA2uni[idx, j]  # differences are at most 0.003 percent!
      decorBDA1uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp1)/sum_wt) 
      decorBDA2uni[idx, j] = 100.0 * (1.0 -numpy.sum(weight * temp2)/sum_wt)

  endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print( "End at %s" % endtime)
  process_end = time.time()
  duration = (process_end - process_start)/3600.0
  print( "pre-plot Total run time: %7.2f hours" % duration)

  x = LHAh
  plt.xlim(-6, 6)
  colours = ['b', 'r', 'b', 'r', 'g','b', 'r', 'b', 'r', 'g']
  colours = ['b', 'r', 'b', 'r', 'g']
  colours = ['b', 'r', 'b', 'r', 'b']
  for i in range(len(dec)):
    plt.plot(x, decorBDA1uni[i,:], colours[i])
    plt.plot(x, decorBDA2uni[i,:], colours[i+1])
  plt.xlabel('Hour Angle (hours)')
  plt.ylabel('Maximum Difference (percent)')
  plt.title('Uniform Wt Theoretical Analysis vs SKA Simulator with Scheme 1 BDA')
  plot_file = 'decor_uniforml_BDA_decorr'
# read in actual data
  print( 'processing data file: ', argv[1])
  dec_values = BDA_utility_functions.getdata(argv[1], argv[2])
  ha =  [-4.5, -3, -1.5, 0, 1, 2, 4]
  ecs = ['10', '-10', '-30', '-50', '-70','10', '-10', '-30', '-50', '-70']
  ecs = ['-70 l', '-70 m']
  plot_parm = ['g^', 'ro']
  plot_parm = ['bo', 'r^', 'bs', 'r^', 'g^','bo', 'ro', 'bs', 'r^', 'g^']
  for i in range(len(dec_values)):
      x = ha
      y = dec_values[i]
      if i < 5:
        plt.plot(x, y,plot_parm[i], label= ecs[i])
      else:
        plt.plot(x, y,plot_parm[i])
# legend(loc=3)
# legend(loc=2)
# legend(loc=4)
# plt.legend(loc=1)
  plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
  plot_file = 'wijnholds_uniform_weighting_double_BDA_plot'
  plot_file = 'Figure_7_plot.png'
  print( 'saving figure ', plot_file)
  plt.savefig(plot_file)
  plt.show()
  plt.clf()
  plt.cla()
  plt.close()

#=============================
# argv[1]  statistics file
# argv[2]  data selection descriptor
if __name__ == "__main__":
  main(sys.argv)
