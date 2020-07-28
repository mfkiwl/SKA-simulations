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

  print(' ')
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  process_start = time.time()
  print("Start at %s" % startime)

  Lx, Ly, Lz = BDA_utility_functions.get_ant_pos()     # get antenna x, y, z ITRF positions
  L = numpy.sqrt(Lx*Lx + Ly*Ly + Lz*Lz)
# Luv = numpy.sqrt(Lx*Lx + Ly*Ly) 
  num_baselines = L.shape[0]

## input parameters for decorrelation calculations
  pi = math.pi
  freq = 700e6                # frequency in MHz
  c = 2.99792e8               # speed of light in m/s
  wavelength = c / freq       # wavelength in m
  dec = numpy.arange(-70,20,20)           # assumed declinations in degrees
  dec = numpy.radians(dec)
 #dec = [-70.0, 10.0]
  LHAh = numpy.arange(-6,6.1,2)   # local hour angle range in hours (2 hr increment)
  LHAh = numpy.arange(-6,6.1,0.1) # local hour angle range in hours (0.1 hr increment)
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

## calculate decorrelation assuming natural weighting
  rotationE = omegaE / wavelength
  pi_t = pi * T
  pi_t_const = 100.0 * (pi_t * pi_t /6)

# calculate decorrelation weights without BDA for uniform weighting
 # weights approximuated based on radial density of (u, v)-points
  threshhold = -1
  threshhold = 50
  w_factor_weight = 0.15
  if single_proc:
    weight_uni = BDA_utility_functions.calculate_basic_uniform_BDA_weight_single(Lx,Ly,Lz,dec,GMST,wavelength, threshhold)
  else:
     weight_uni,weight_uni1 = BDA_utility_functions.calculate_basic_uniform_BDA_weight(Lx,Ly,Lz,dec,GMST,wavelength,threshhold,w_factor_weight)

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
      weight = weight_uni[idx,j,:]
      sum_wt = numpy.sum(weight)
      inv_sum = pi_t_const / sum_wt
    # calculate decorrelation
#     temp1 = numpy.square(dudt * l1 + dvdt * m1 );
#     temp2 = numpy.square(dudt * l2 + dvdt * m2 );
      temp1 = numpy.square(dudt * l1 + dvdt * m1 + dwdt *n1);
      temp2 = numpy.square(dudt * l2 + dvdt * m2 + dwdt *n2);
      decorBDA1uni[idx, j] = inv_sum * numpy.sum(weight * temp1) 
      decorBDA2uni[idx, j] = inv_sum * numpy.sum(weight * temp2)

  endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("End at %s" % endtime)
  process_end = time.time()
  duration = (process_end - process_start)/3600.0
  print("pre-plot Total run time: %7.2f hours" % duration)

  x = LHAh
  plt.xlim(-6, 6)
  colours = ['g', 'r', 'b', 'r', 'b']
  print('declinations', numpy.degrees(dec))
  for i in range(len(dec)):
    if i == 0:
      print('dec colour', numpy.degrees(dec[i]), colours[i])
      print('l, decor_l', decorBDA1uni[i,:])
      plt.plot(x, decorBDA1uni[i,:], colours[i])
      print('m, decor_m', decorBDA2uni[i,:])
      plt.plot(x, decorBDA2uni[i,:], colours[i])
    if i == len(dec)-1:
      print('dec colour', numpy.degrees(dec[i]), colours[i])
      print('l, decor_l', decorBDA1uni[i,:])
      plt.plot(x, decorBDA1uni[i,:], colours[i])

  plt.xlabel('Hour Angle (hours)')
  plt.ylabel('Maximum Difference (percent)')
  plt.title('Uniform Wt Theoretical Analysis vs SKA Simulator with 0.14 sec Integration')
# read in actual data
  try:
    print('processing data file: ', argv[1])
    dec_values = BDA_utility_functions.getdata(argv[1], argv[2])
    good_dec = True
  except:
    good_dec = False
# ha =  [-4.5, -3, -1.5, 0, 1, 2, 4]
# ecs = ['10', '-10', '-30', '-50', '-70','10', '-10', '-30', '-50', '-70']
# ecs = ['10', '-70', '-70']
# plot_parm = ['bo', 'g^', 'g^']

  ha =  [-4.5, -3, -1.5, 0, 1, 2, 4]
  ecs = ['10', '-10', '-30', '-50', '-70','10', '-10', '-30', '-50', '-70']
  plot_parm = ['bo', 'ro', 'bs', 'r^', 'g^','bo', 'ro', 'bs', 'r^', 'g^']

  if good_dec:
    print('len(dec_values)', len(dec_values))
    print(dec_values)
    for i in range(len(dec_values)):
        x = ha
        y = dec_values[i] 
        try:
          if i < len(dec_values) -1:
            plt.plot(x, y,plot_parm[i], label= ecs[i])
          else:
            plt.plot(x, y,plot_parm[i])
        except:
          continue
# plt.legend(loc=1)
  plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
  plot_file = 'sjw_uniform_weighting_007_014sec_avg_plot_100_with_w'
  plot_file = 'Figure_3b_plot.png'
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
