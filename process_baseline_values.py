import os
import time
import sys
import numpy
import pyrap.tables
import math
from multiprocessing import Process, Queue

try:
   import cPickle as pickle
except:
   import pickle

def baseline_worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)


def process_baseline(num_avg, array_to_process, UVWs, key):
  num_points = 0
  max_i = array_to_process.shape[0]
  for i in range(0,max_i,num_avg):
    num_points = num_points + 1
    sum = 0.0
    sum_uvw = 0.0
    num_sum = 0
    for j in range(i, i+num_avg):
      if j < max_i:
        sum = sum + array_to_process[j,:,:]
        sum_uvw = sum_uvw + UVWs[j,:]
        num_sum = num_sum + 1
    avg = sum / num_sum
    avg_uvw = sum_uvw / num_sum
    for j in range(i, i+num_avg):
      try:
        if j < max_i:
          array_to_process[j,:,:] =  avg[:,:]
          UVWs[j,:] = avg_uvw[:]
      except:
        pass
# print 'processed array ', key, num_avg, num_points
# print 'shapes ', avg.shape, avg_uvw.shape
  return max_i, num_points, key, array_to_process, UVWs

def analyse_baselines(data_ms,data_values,data_flag,salvini_flag):
  data_flag = int(data_flag)
  salvini_flag = int(salvini_flag)

# read in pickle file with information about baseline lengths
  infile = open(data_values, 'rb')
  baseline_data = pickle.load(infile, encoding='bytes')
  print('number of keys ', len(baseline_data))
  min_value = 1.0e6
  max_value = 0.0
  keys = baseline_data.keys()
  for key in keys:
    if baseline_data[key][0] < min_value:
      min_value =  baseline_data[key][0]
      min_key = key
    if baseline_data[key][0] > max_value:
      max_value =  baseline_data[key][0]
      max_key = key

  print('min and max keys ', min_key, max_key)
  print('min and max baselines ', min_value, max_value)
# baselines 160, 80, 40, 20, 10, 5 km zones
# avg         0   2   4   8  16  32
  baseline_zones = [30,15,7.5,3.75,1.875,0.9375,0.46875,0.234375]
  baseline_avg = [0,6,12,24,48,96,192,384,768]

# tried the following - result gives greater residuals than prev
# and not much more compression
# baseline_avg = [0,7,15,31,62,112,225,450]

# checking on this (averaging method 1)
  baseline_zones = [80,40,20,10,5,2.5,1.25,0.625,0.3125]
  baseline_avg = [1,2,4,8,16,32,64,128,256,512]

  baseline_zones = [80, 40, 30, 20, 15, 10, 7.5, 5, 3.75, 2.5, 1.875, 1.25, 0.9375, 0.625, 0.5625, 0.375, 0.28125]
  baseline_avg = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]

  L = len(baseline_zones)

# open MS
  t = pyrap.tables.table(data_ms,readonly=False)
  if data_flag == 2:
    print(' ++++ averaging CORRECTED_DATA')
    data = t.getcol('CORRECTED_DATA')
  if data_flag == 1:
    print(' ++++ averaging MODEL_DATA')
    data = t.getcol('MODEL_DATA')
  if data_flag == 0:
    print(' ++++ averaging DATA')
    data = t.getcol('DATA')

# store output in MODEL_DATA column
# data_flag = 1  

  uvws = t.getcol('UVW')
  print('uvws shape ', uvws.shape)

  print('salvini flag is ', salvini_flag)
  process_num = 0
  initial_points = 0
  final_points = 0
  num_long_baselines = 0
  eps = numpy.spacing(1)
  maxbmult = 125
  maxbmult = 500
  maxbmult = 63
  maxbmult = 32
  maxbmult = 250
  bmax = max_value / 1000.0   # maximum baseline in km
  TASKS = []
  for key in keys:
    initial_points = initial_points + len(baseline_data[key]) - 1
    baseline_length = baseline_data[key][0] / 1000.0 # baseline length in km
    if salvini_flag == 0:
      num_avg = baseline_avg[L]
      for i in range(L):
        if baseline_length  > baseline_zones[i]:
          num_avg = baseline_avg[i]
          break
    else:
# Following line give the Stef Salvini algorithm for determining the number of points to average
      num_avg = min(maxbmult,int(numpy.fix((bmax+100*eps)/baseline_length)));
      print('key salvini num_avg ', key, num_avg)
    if num_avg > 1:
      array_to_process = data[baseline_data[key][1:]]
      UVWs = uvws[baseline_data[key][1:]]
      TASKS.append((process_baseline, (num_avg, array_to_process, UVWs, key)))
    else:
      num_long_baselines = num_long_baselines + 1
      process_num = process_num + 1
      num_points = len(baseline_data[key]) - 1
      final_points = final_points + num_points
 
  print('number of long baselines ', num_long_baselines)
  print('number of averaging tasks ', len(TASKS))

  # Create queues
  task_queue = Queue()
  done_queue = Queue()

  # Submit tasks
  for task in TASKS:
     task_queue.put(task)

  # Start worker processes
  num_processors = 1
  if num_processors <= 2:
    try:
      import multiprocessing
      processors =  multiprocessing.cpu_count()
      if processors > num_processors:
        num_processors = processors
        print('*** setting number of processors to',num_processors)
    except:
      pass
  num_processors = 1

  for i in range(num_processors):
     Process(target=baseline_worker, args=(task_queue, done_queue)).start()

  # Get and print results
  for i in range(len(TASKS)):
    init_points, num_points, key, modified_array, modified_uvws  = done_queue.get(timeout=300)
#   print 'for baseline ', key, 'initial and final number of points ', init_points, num_points, float(init_points) / num_points

    data[baseline_data[key][1:]] = modified_array
    uvws[baseline_data[key][1:]] = modified_uvws
    final_points = final_points + num_points
    process_num = process_num + 1
    if process_num % 1000 == 0:
      print('processed baseline number ', process_num)

  # Tell child processes to stop
  for i in range(num_processors):
    task_queue.put('STOP')

  if data_flag == 2:
    print( ' ++++ storing CORRECTED_DATA')
    t.putcol('CORRECTED_DATA', data)
  if data_flag == 1:
    print( ' ++++ storing MODEl_DATA')
    t.putcol('MODEL_DATA', data)
  if data_flag == 0:
    print( ' ++++ storing DATA')
    t.putcol('DATA', data)
  print (' ++++ storing averaged UVWs')
  t.putcol('UVW', uvws)
  t.flush()
  t.close()         
  print('finished processing')
  print('initial and final points ', initial_points, final_points)
  if salvini_flag:
    print('Salvini percentage decrease ', ((initial_points - final_points) * 100.0) / initial_points)
  else:
    print('Willis percentage decrease ', ((initial_points - final_points) * 100.0) / initial_points)

def main( argv ):
# argv[1] = name of MS 
# argv[2] = name of pickle file with baseline info
# argv[3] = flag to get MODEL_DATA or DATA column
# argv[4] = flag to use Salvini or Willis algorithm
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("Start at %s" % startime)

  analyse_baselines(argv[1], argv[2], argv[3], argv[4])

  endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("End at %s" % endtime)
  return
#=============================
if __name__ == "__main__":
  main(sys.argv)

