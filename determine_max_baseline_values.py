#!/usr/bin/env python

# script to analyze Measurement Set baseline parameters and store them in a pickle file

import os
import time
import sys
import numpy 
import pyrap.tables
import math 

try:
   import cPickle as pickle
except:
   import pickle

# function that does the 'real' work
def analyse_baselines(data_ms, use_just_uv):

# open MS
  t = pyrap.tables.table(data_ms,readonly=True)

  # get U,V,W s from MS
  uvws = t.getcol('UVW')
  ant1 = t.getcol('ANTENNA1')
  ant2 = t.getcol('ANTENNA2')
  times = t.getcol('TIME')
  scan_number = t.getcol('SCAN_NUMBER')
  spwid = t.getcol('DATA_DESC_ID')
  u = uvws[:,0] 
  v = uvws[:,1]
  w = uvws[:,2]
  if use_just_uv == 1:
    baseline = numpy.sqrt(u*u + v*v)
  else:
    baseline = numpy.sqrt(u*u + v*v + w*w)
   
  print('baseline shape', baseline.shape)
  time = 0.0
  scan = 0
  baseline_group = -1
  baseline_groups = {}
  scan_starts = []
  scan_starts.append(scan)
  for i in range(scan_number.shape[0]):
    if (scan_number[i] - scan) > 0:
      scan = scan_number[i]
      scan_starts.append(scan)
      print('new scan at time ', times[i])
      max_baseline = 0.0
      baseline_group = baseline_group + 1
      print('new baseline group is ',  baseline_group)
    ants = (ant1[i], ant2[i])
    if baseline[i] > max_baseline:
      max_baseline = baseline[i]
      max_ants= ants
    dict_key = (spwid[i],baseline_group, ants)
    if dict_key in baseline_groups:
      baseline_groups[dict_key].append(i)
      if baseline[i] > baseline_groups[dict_key][0]:
        baseline_groups[dict_key][0] =  baseline[i] 
    else:
      baseline_groups[dict_key] =  [baseline[i], i]
  scans_end =  scan_number.shape[0]
  scan_starts.append(scans_end)
  print('\nfinished baseline analysis!!')
  # save output in pickle format
  print('longest baseline ', max_baseline, max_ants)
  keys =  baseline_groups.keys()
  if use_just_uv == 1:
    out_pickle = data_ms + '_max_baselines_uv_values'
  else:
    out_pickle = data_ms + '_max_baselines_values'
  outputx = open(out_pickle,'wb')
  pickle.dump(baseline_groups, outputx)
  outputx.close()

  out_pickle = data_ms + '_scan_positions'
  outputx = open(out_pickle,'wb')
  pickle.dump(scan_starts, outputx)
  outputx.close()

  t.close()         


def main( argv ):
  print('starting baseline analysis')
# argv[1] = name of measurement set
# argv[2] = use just uv or u,v,w
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("Start at %s" % startime)
  try:
    use_just_uv = int(argv[2])
  except:
    use_just_uv = 0
  print('data selection for uv = ', use_just_uv)

  analyse_baselines(argv[1],use_just_uv)
  endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("End at %s" % endtime)
  return

#=============================
if __name__ == "__main__":
  main(sys.argv)
