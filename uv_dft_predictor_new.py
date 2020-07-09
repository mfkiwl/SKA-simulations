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

###########################################
def lm_to_radec_deg(ra0,dec0,l,m):
# using SIN projection

    sind0=math.sin(dec0)
    cosd0=math.cos(dec0)
    n = math.sqrt(1-l**2-m**2)
    dec = math.asin(m * cosd0 + n * sind0)

    bottom = n * cosd0 - m * sind0
    ra = ra0 + math.atan2(l, bottom)

    dec_deg = math.degrees(dec)
    ra_deg = math.degrees(ra)
    return (ra_deg,dec_deg)

###########################################
def deg_lm_to_radec(ra0,dec0,l,m):
# assume inputs are in degrees
# using SIN projection

    el = math.radians(l)
    em = math.radians(m)
    dec_0 = math.radians(dec0)
    ra_0 = math.radians(ra0)
    sind0=math.sin(dec_0)
    cosd0=math.cos(dec_0)
    en = math.sqrt(1-el**2-em**2)
    dec_deg = math.degrees(math.asin(em * cosd0 + en * sind0))

    bottom = en * cosd0 - em * sind0
    ra_dec = math.degrees(ra_0 + math.atan2(el, bottom))
#   returned values are in degrees
    return (ra_deg,dec_deg)

###########################################
def radec_to_lmn(ra0,dec0,ra,dec):
    l=math.cos(dec)*math.sin(ra-ra0)
    sind0=math.sin(dec0)
    m=math.sin(dec) * math.cos(dec0) - math.cos(dec) * math.sin(dec0) * math.cos(ra - ra0)  # from pandey;  gives same results for casa and cyga
    n = math.sqrt(1-l**2-m**2)
    return (l,m,n)

###########################################
def deg_radec_to_lmn(ra0,dec0,ra,dec):
    ra_0 = math.radians(ra0)
    dec_0 = math.radians(dec0)
    ra_r = math.radians(ra)
    dec_r = math.radians(dec)
    l=math.cos(dec_r)*math.sin(ra_r-ra_0)
    sind0=math.sin(dec_0)
    m=math.sin(dec_r) * math.cos(dec_0) - math.cos(dec_r) * math.sin(dec_0) * math.cos(ra_r - ra_0)  
    n = math.sqrt(1-l**2-m**2)
    l_deg = math.degrees(l)
    m_deg = math.degrees(m)
    n_deg = math.degrees(n)
    return (l_deg,m_deg,n_deg)


# compute visibilities by dft
def compute_uv_values(uvws, source_list,uv_scale_factor,num_avg): 
    returned_data = []
    returned_data.append(uv_scale_factor.shape[0])
    for i in range(len(uv_scale_factor)):
      returned_data.append(i)   # equates to freq
      uvw_group = uvws * uv_scale_factor[i]

      u = uvw_group[:,0]
      v = uvw_group[:,1] 
      w = uvw_group[:,2]

      if num_avg > 1:
        num_uvws = u.shape[0]
        x , step = numpy.linspace(0, num_uvws, num = num_uvws * num_avg, endpoint=False,retstep=True)
        xp, step1 = numpy.linspace(0, num_uvws, num = num_uvws, endpoint=False,retstep=True)
        u_interp = numpy.interp(x, xp, u)
        v_interp = numpy.interp(x, xp, v)
        w_interp = numpy.interp(x, xp, w)
      else:
        u_interp = u
        v_interp = v
        w_interp = w

      for j in range(len(source_list)):
        source_parm = source_list[j]
        flux = source_parm[2]
        l = source_parm[0]
        m = source_parm[1]
        n = math.sqrt(1 - l*l - m*m)
        phase = 2 * numpy.pi * (u_interp * l + v_interp * m + w_interp * (n - 1))
        if j == 0:
          uv_data_vals = flux * numpy.exp(phase* 1j)
        else:
          uv_data_vals = uv_data_vals + flux * numpy.exp(phase* 1j)
#     print('raw averaged uv_data shape, vals', uv_data_vals.shape, uv_data_vals)
      if num_avg > 1:
        offset = int(math.ceil(num_avg / 2.0))
        avg_list = []
        avg_list.append(0.5 * uv_data_vals[0] + 0.5 * numpy.average(uv_data_vals[1:offset]))
        for j in range (num_avg-offset+1, uv_data_vals.shape[0]-num_avg, num_avg):
          avg_list.append(numpy.average(uv_data_vals[j:j+num_avg]))
        uv_data_vals = numpy.array(avg_list)
        abs_vals = numpy.absolute(uv_data_vals)
#       print 'min and max ', numpy.amin(abs_vals), numpy.amax(abs_vals)
#       uv_data_vals = uv_data_vals / abs_vals
        
      returned_data.append(uv_data_vals)
#     returned_data.append(uv_data_vals[0:num_uvws])
#   print'returning from compute_uv_values'
    return returned_data

def process_baseline(array_to_process, UVWs, key, source_list, uv_scale_factor, num_avg):
  result = compute_uv_values(UVWs,source_list,uv_scale_factor, num_avg)
  num_freq = result[0]
  index = 0
  for k in range(num_freq):
    index = index + 1
    uv_freq = result[index]
    if uv_freq != k:
      print('discrepancy in data !!!!')
    index = index + 1
    uv_data_value = result[index]
    array_to_process[:,k,0] = uv_data_value
    array_to_process[:,k,3] = uv_data_value
    array_to_process[:,k,1] = 0+0j
    array_to_process[:,k,2] = 0+0j
  return key, array_to_process

def analyse_baselines(data_ms,data_values,data_flag, source_offset, num_sub,dir_offset):

# argv[1] = name of MS 
# argv[2] = name of pickle file with baseline info
# argv[3] = flag to get MODEL_DATA or DATA column
# argv[4] = source offset in arcmin
# argv[5] = number of sub averaging intervals
# argv[6] = offset is in l (0) or m (1)

  os.system('date')
  process_start = time.time()
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("uv_dft_predictor  Start at %s" % startime)

# open MS
  t = pyrap.tables.table(data_ms,readonly=False)
  print('opened MS')

  direction = pyrap.tables.table(t.getkeyword('FIELD'))
  ref_direction = direction.getcol('REFERENCE_DIR')
  ra0 =  ref_direction[0][0][0]
  dec0 =  ref_direction[0][0][1]
  print('field centre (hours, deg)', (24.0/360.0) * math.degrees(ra0), math.degrees(dec0))

# get frequencies from MS
  spw = pyrap.tables.table(t.getkeyword('SPECTRAL_WINDOW'))
  frequency = spw.getcol('CHAN_FREQ')[0]
  print('*** frequency is', frequency)
  speed_of_light = 299792458     # metres / sec
  uv_scale_factor = frequency / speed_of_light  # wavelength (in metres)^-1



  # create some dummy sources
# single source at specified offset
  diag = 1.0
  source_list = []
  source_offset = float(source_offset)
  dist = math.radians(source_offset / 60.0)
  source_m = dist 
  source_l = 0.0 
  source_flux = 1.0
  source_parm = (source_l, source_m, source_flux)
  print('field centre (radians) ', ra0, dec0)
  print('l, m (radians) ', source_l, source_m)
  ra_loc, dec_loc = lm_to_radec_deg(ra0,dec0,source_l,source_m)
  print('source position ', ra_loc, dec_loc)
  print('source_parm', source_parm)
  source_list.append(source_parm)

  dist = math.radians(source_offset / 60.0 + 9.0 / 3600)
  source_m =  dist
  source_l = math.radians(9.0 / 3600.0)
  source_flux = 1.0 / 1000000.0
  source_parm = (source_l, source_m, source_flux)
  ra_loc, dec_loc = lm_to_radec_deg(ra0,dec0,source_l,source_m)
  print('source position ', ra_loc, dec_loc)
  source_list.append(source_parm)

  source_l = math.radians(-9.0 / 3600.0)
  source_parm = (source_l, source_m, source_flux)
  ra_loc, dec_loc = lm_to_radec_deg(ra0,dec0,source_l,source_m)
  print('source position ', ra_loc, dec_loc)
  source_list.append(source_parm)

  dist = math.radians(source_offset / 60.0 - 9.0 / 3600)
  source_m = dist
  source_l = math.radians(9.0 / 3600.0)
  source_flux = 1.0 / 100000.0
  source_parm = (source_l, source_m, source_flux)
  ra_loc, dec_loc = lm_to_radec_deg(ra0,dec0,source_l,source_m)
  print('source position ', ra_loc, dec_loc)
  source_list.append(source_parm)
  
  source_l = math.radians(-9.0 / 3600.0)
  source_parm = (source_l, source_m, source_flux)
  ra_loc, dec_loc = lm_to_radec_deg(ra0,dec0,source_l,source_m)
  print('source position ', ra_loc, dec_loc)
  source_list.append(source_parm)

# single source at specified offset
  source_list = []
  source_offset = float(source_offset)
  dist = math.radians(source_offset / 60.0)
  print('direction offset ', int(dir_offset))
  if int(dir_offset) > 0:
    source_m = dist 
    source_l = 0.0 
    print('source offset in m = ', source_m )
  else:
    source_l = dist 
    source_m = 0.0 
    print('source offset in l = ', source_l )
  source_flux = 1.0
  source_parm = (source_l, source_m, source_flux)
  print('field centre (radians) ', ra0, dec0)
  print('l, m (radians) ', source_l, source_m)
  ra_loc, dec_loc = lm_to_radec_deg(ra0,dec0,source_l,source_m)
  print('****** source position ', ra_loc, dec_loc)
  print('source_parm', source_parm)
  source_list.append(source_parm)

  data_flag = int(data_flag)
  num_avg = int(num_sub)
  test_val = int(math.ceil(num_avg / 2.0))
  if test_val == num_avg / 2:
    print('number for interpolation must be odd number so adding one to ', num_avg)
    num_avg = num_avg +1
    print('number to interpolate set to ', num_avg)

# read in pickle file with information about baseline lengths
  openfile = open(data_values, 'rb')
  baseline_data = pickle.load(openfile,encoding='bytes')
  print('number of keys ', len(baseline_data))
  keys = baseline_data.keys()
  openfile.close()

# open MS - can get data from any column - we jsut use CORRECTED_DATA for convenience
  data = t.getcol('CORRECTED_DATA')
  uvws = t.getcol('UVW')
  print('uvws shape ', uvws.shape)

  TASKS = []
  for key in keys:
      array_to_process = data[baseline_data[key][1:]]
      UVWs = uvws[baseline_data[key][1:]]
      TASKS.append((process_baseline, (array_to_process, UVWs,key, source_list, uv_scale_factor, num_avg)))

  print('number of processing tasks ', len(TASKS))

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

# num_processors = 1
  for i in range(num_processors):
     Process(target=baseline_worker, args=(task_queue, done_queue)).start()

  process_num = 0
  # Get and print results
  for i in range(len(TASKS)):
    key, modified_array = done_queue.get(timeout=300)
    data[baseline_data[key][1:]] = modified_array
    process_num = process_num + 1
    if process_num % 1000 == 0:
      print('computed visibilities for baseline number ', process_num)

  # Tell child processes to stop
  for i in range(num_processors):
    task_queue.put('STOP')

  if data_flag == 3:
    print(' ++++ storing CORRECTED_DATA and MODEL_DATA')
    t.putcol('CORRECTED_DATA', data)
    t.putcol('MODEL_DATA', data)
  if data_flag == 2:
    print(' ++++ storing CORRECTED_DATA')
    t.putcol('CORRECTED_DATA', data)
  if data_flag == 1:
    print(' ++++ storing MODEl_DATA')
    t.putcol('MODEL_DATA', data)
  if data_flag == 0:
    print(' ++++ storing DATA')
    t.putcol('DATA', data)
  t.flush()
  t.close()         
  print('finished processing')

def main( argv ):
# argv[1] = name of MS 
# argv[2] = name of pickle file with baseline info
# argv[3] = flag telling system which column to store data in
# argv[4] = source offset in arcmin
# argv[5] = number of sub averaging intervals
# argv[6] = offset is in l (0) or m (1)
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  process_start = time.time()
  print("Start at %s" % startime)

  analyse_baselines(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6])

  endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("End at %s" % endtime)
  process_end = time.time()
  duration = (process_end - process_start)/3600.0
  print("uv_dft_predictor Total run time: %7.2f hours" % duration)

  return
#=============================
if __name__ == "__main__":
  main(sys.argv)


