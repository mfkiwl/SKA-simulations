import numpy
import math
from multiprocessing import Process, Queue

def baseline_worker(input, output):
  for func, args in iter(input.get, 'STOP'):
    result = func(*args)
    output.put(result)

# Convert latitude (radians), longitude (radians) and elevation (metres) to ITRF XYZ
def WGS84ToITRF(lat, lon, h): # WGS-84 to ITRF

# Earth oblateness constants
  sm_a = 6378137.0
  invf = 298.257223563
  f = 1.0 / invf

  SINK = math.sin(lat)
  COSK = math.cos(lat)
  e2 = 2.0 * f - f * f
  v = sm_a / math.sqrt(1.0 - e2 * SINK * SINK)
  x = (v + h) * COSK * math.cos(lon)
  y = (v + h) * COSK * math.sin(lon)
  z = ((1 - e2) * v + h) * SINK
  return x, y, z

# get actual SKA test observation data
def getdata(filename, fits_text):
    print('parameters', filename, fits_text)
    text = open(filename, 'r').readlines()
    L = len(text)
    j = 0
    dec_values = []
    values = []
    counter = 0
# data = (fitsimage,sigma,mindata,minlocation[1][0],minlocation[0][0],maxdata,maxlocation[1][0],maxlocation[0][0],((maxdata-mindata)/sigma))

    for i in range(L):
#      try:
          location_fits = text[i].find(fits_text)
          if location_fits > -1:
#             print('processing ', text[i])
              info = text[i].split()
              print('info', info)
              if info[0] == '#':
                continue
              parm = 100.0 * (1.0 - float(info[5]))
#             parm = (1.0 - float(info[5]))
#             parm = float(info[5])
              counter = counter + 1
              values.append(parm)
              j = j + 1
              if j > 6:
                 dec_values.append(values)
                 values = []
                 j = 0
#      except:
#         pass
#   print 'dec_values', dec_values
    return dec_values 
              

# get SKA telescope positions from text file
def get_ant_pos(convert_to_itrf=True):
  xx=[]
  yy=[]
  zz=[]
  print('loading antenna positions')
  if convert_to_itrf:
    SKA1_filename = 'ska_pos.csv'
  else:
    SKA1_filename = 'vla_a.csv'
  text = open(SKA1_filename, 'r').readlines()
  L = len(text)
  for i in range(L):
    try:
       info = text[i].split()
#      print 'info ', info
       if convert_to_itrf:
           lon = math.radians(float(info[0]))
           lat = math.radians(float(info[1]))
           el = float(info[2])
           x, y, z = WGS84ToITRF(lat, lon, el)
       else:
           x = math.radians(float(info[0]))
           y = math.radians(float(info[1]))
           z = math.radians(float(info[2]))
       xx.append(x)
       yy.append(y)
       zz.append(z)
    except:
      pass
  
  xx = numpy.array(xx)
  yy = numpy.array(yy)
  zz = numpy.array(zz)

  num_ants = len(xx)
  print('number of antennas ', num_ants)
  num_baselines = int(0.5 *(num_ants * (num_ants-1)))
  baselines = numpy.zeros((num_baselines, 3), numpy.float64);
  k = 0
  for i in range(num_ants):
    for j in range(i+1,num_ants):
      baselines[k,0] = xx[i] - xx[j]
      baselines[k,1] = yy[i] - yy[j]
      baselines[k,2] = zz[i] - zz[j]
      k = k + 1
  Lx =baselines[:,0]
  Ly =baselines[:,1]
  Lz =baselines[:,2]

  return Lx, Ly, Lz

# compute u,v,w positions and rate of change
def uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin):

  u = Lx * lha_sin + Ly * lha_cos
  v = (-Lx * lha_cos + Ly * lha_sin) * dec_sin + Lz * dec_cos
  u_max = numpy.max(u)
  v_max = numpy.max(v)
# print 'u_max v_max', u_max, v_max
# note adding the Lz termm above causes some sort of round-off error in the 
# calculate_basic_uniform_BDA_weight function
# v = (-Lx * lha_cos + Ly * lha_sin) * dec_sin 
  w = (Lx * lha_cos - Ly * lha_sin) * dec_cos  + Lz * dec_sin

  dudt = (Lx * lha_cos - Ly * lha_sin) * rotationE
  dvdt = (Lx * lha_sin + Ly * lha_cos) * dec_sin * rotationE;
  dwdt = (Lx * lha_sin + Ly * lha_cos) * dec_cos * rotationE * (-1.0)

  return u,v,w,dudt,dvdt,dwdt, u_max, v_max

# calculate weights for uniform gridding equations
# parallel processing version
def calculate_basic_uniform_BDA_weight(Lx,Ly,Lz,dec,GMST,wavelength,threshhold=-1, w_factor_weight=0.15):

  omegaE = 2 * math.pi / (23 * 3600 + 56 * 60)  # angular velocity of Earth rotation
                                         # in  rad/s
  rotationE = omegaE / wavelength

  weight_uni = numpy.zeros((len(dec), len(GMST),len(Lx)), numpy.float64)
  weight_uni1 = numpy.zeros((len(dec), len(GMST),len(Lx)), numpy.float64)

# calculate decorrelation without BDA for uniform weighting
  TASKS = []
  for idx in range(len(dec)):
    dec_sin =  math.sin(dec[idx])
    dec_cos =  math.cos(dec[idx])
    for j in range(len(GMST)):
      lha_cos = math.cos(GMST[j])
      lha_sin = math.sin(GMST[j])
      TASKS.append((process_baseline,(idx,j,threshhold,w_factor_weight,Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)))

  print('number of tasks for parallel processings', len(TASKS))
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
  # Get the results from parallel processing
  for i in range(len(TASKS)):
    idx,j, wt, wt1 = done_queue.get(timeout=300)
    weight_uni[idx,j,:] = wt
    weight_uni1[idx,j,:] = wt1
    process_num = process_num + 1
    if process_num % (num_processors * 2) == 0:
      print('****************** caught data for process ', process_num)

  # Tell child processes to stop
  for i in range(num_processors):
    task_queue.put('STOP')

  return weight_uni, weight_uni1

def process_baseline(idx,j,threshhold,w_factor_weight,Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin):

  u,v,w,dudt,dvdt,dwdt,u_max,v_max = uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)
  len_u = u.shape[0]
  weight = numpy.zeros((len_u,), numpy.float64)
  weight1 = numpy.zeros((len_u,), numpy.float64)
  N = 1.0
  for k in range(len(u)):
    weight1[k] = numpy.sqrt(dudt[k] * dudt[k] + dvdt[k] * dvdt[k] + w_factor_weight * (dwdt[k] * dwdt[k])) 
    if threshhold > 0:
      u_sq = numpy.square(u - u[k]) # note Lz* dec_cos term should cancel here 
      v_sq = numpy.square(v - v[k])
      w_sq = numpy.square(w - w[k])
      uvdist = numpy.sqrt(u_sq + v_sq + w_sq)
      N = numpy.sum(uvdist < threshhold) * 1.0
    div = 1.0 / N
    weight[k] = div * weight1[k]
#   if N > 1:
#     print('dec_index,GMST_index,weight1[k],div,k,N', idx,j,weight1[k],div,k, N)
  return idx,j,weight, weight1

# calculate weights for uniform gridding equations
# single processor version
def calculate_basic_uniform_BDA_weight_single(Lx,Ly,Lz,dec,GMST,wavelength):

  omegaE = 2 * math.pi / (23 * 3600 + 56 * 60)  # angular velocity of Earth rotation
                                         # in  rad/s
  rotationE = omegaE / wavelength

# calculate decorrelation without BDA for uniform weighting
  weight_uni = numpy.zeros((len(dec), len(GMST),len(Lx)), numpy.float64);
  threshhold = 150             # distance threshold to determine density 
  threshhold = 300             # distance threshold to determine density 
  threshhold = 200             # distance threshold to determine density (seems to give 'best' results) 
  for idx in range(len(dec)):
    print('computing weights for dec', math.degrees(dec[idx]), ' ...')
    dec_sin =  math.sin(dec[idx])
    dec_cos =  math.cos(dec[idx])
    for j in range(len(GMST)):
      lha_cos = math.cos(GMST[j])
      lha_sin = math.sin(GMST[j])
      u,v,w,dudt,dvdt,dwdt,u_max,v_max = uv_plane_parms(Lx,Ly,Lz,rotationE,dec_sin,dec_cos,lha_cos,lha_sin)
      for k in range(len(u)):
#       u_sq = numpy.square(u - u[k]) # note Lz* dec_cos term should cancel here 
#       v_sq = numpy.square(v - v[k])
#       uvdist = numpy.sqrt(u_sq + v_sq)
#       N = numpy.sum(uvdist < threshhold)
#       weight_uni[idx,j,k] = (1.0/N) * numpy.sqrt(dudt[k] * dudt[k] + dvdt[k] * dvdt[k]) 
        weight_uni[idx,j,k] = u_max * v_max * numpy.sqrt(dudt[k] * dudt[k] + dvdt[k] * dvdt[k]) 
#       weight_uni[idx,j,k] = u_max * v_max * numpy.sqrt(dudt[k] * dudt[k] + dvdt[k] * dvdt[k] + dwdt[k] * dwdt[k]) 
    print('computed weights for dec', math.degrees(dec[idx]), ' ...')

  return weight_uni
