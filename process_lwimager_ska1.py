# script to image the calibrated visibility data 

import os
import time
import sys

def lwimagerMS(msName,imageName,stokes_parm,index):

# delete any stuff left from previous run
   fits_image = imageName + "_" + stokes_parm + "_"+str(i)+".fits"
   os.system("/bin/rm -rf " + fits_image)
   remove_file = msName + "-" + stokes_parm + "-*"
   command = '/bin/rm -rf ' + remove_file
   os.system(command)
   os.system("/bin/rm -rf TempLattice*")

# build command
   cc = 'lwimager '
   cc += ' ms='+msName +'.MS'
   cc += ' mode=channel'
   cc += ' nchan=1'
   cc += ' chanstep=1'
   cc += ' chanstart='+str(index) 
   cc += ' img_nchan=1'
   cc += ' img_chanstep=1'
   cc += ' img_chanstart='+str(index) 
   if fits_image.find('natural') > 0:
     cc += ' weight=natural'
   else:
     cc += ' weight=uniform'
   cc += ' data=CORRECTED_DATA'
   cc += ' stokes='+stokes_parm
   cc += ' cellsize=0.25arcsec'
#  cc += ' cellsize=0.10arcsec'
   cc += ' npix=1024'
   cc += ' fits=' + fits_image

#  cc += ' niter=1'             # lwimager seems to insist on at least one clean component
   cc += ' niter=10'
   cc += ' gain=0.1'
   cc += ' operation=csclean'
   cc += ' threshold=0mJy'


   if fits_image.find('2m') > 0:
     phase_str = " phasecenter='j2000, 0.00416666666667deg, -7.9995936206deg' "
# M offset 2 deg
     if fits_image.find('_70') > 0:
       phase_str = " phasecenter='j2000, 0.00416666666667deg, -67.9995936206deg' "
     if fits_image.find('_50') > 0:
       phase_str = " phasecenter='j2000, 0.00416666666667deg, -47.9995936206deg' "
     if fits_image.find('_30') > 0:
       phase_str = " phasecenter='j2000, 0.00416666666667deg, -27.9995936206deg' "
     if fits_image.find('_p10') > 0:
       phase_str = " phasecenter='j2000, 0.00416666666667deg, 12.0004063794deg' "
      

# L offset 2 deg
   if fits_image.find('2l') > 0:
     phase_str = " phasecenter='j2000, 2.03540688472deg, -9.99384320951deg' "
     if fits_image.find('_70') > 0:
       phase_str = " phasecenter='j2000, 5.8351270824deg, -69.9042853228deg' "
     if fits_image.find('_50') > 0:
       phase_str = " phasecenter='j2000, 3.11445288709deg, -49.9584052617deg' "
     if fits_image.find('_30') > 0:
       phase_str = " phasecenter='j2000, 2.31372429234deg, -29.9798425777deg' "
     if fits_image.find('_p10') > 0:
       phase_str = " phasecenter='j2000, 2.03540688472deg, 9.99384320951deg' "

#  phase_str = " phasecenter='j2000, 0.00416666666667deg, -9.2302216307deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, -63.217503809deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, -56.2689405328deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, -41.658257056deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, 18.341742944deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, -21.658257056deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, 38.341742944deg ' "
#  phase_str = " phasecenter='j2000, 40.0053080796deg, -42.3936091359deg ' "
#  phase_str = " phasecenter='j2000, 31.9199632886deg, -26.1080834043deg ' "
#  phase_str = " phasecenter='j2000, 28.7139316947deg, -8.79115455438deg ' " 
#  phase_str = " phasecenter='j2000, 28.7139316947deg, 8.79115455438deg ' "
#  phase_str = " phasecenter='j2000, 57.6256610444deg, -55.7972183277deg ' "
#  phase_str = " phasecenter='j2000, 0.00416666666667deg, -41.658257056deg ' "

   cc += phase_str

   print('using command ', cc)
   os.system(cc)

#=============================
if __name__ == "__main__":
# argv[1] = name of measurement set
# argv[2] = name of output fits file
# argv[3] = hour angle number
# argv[4] = number of channels
  os.system('date')
  process_start = time.time()
  startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("process_lwimager_ska Start at %s" % startime)

  msName=sys.argv[1]
  imageName=sys.argv[2] + '_' + sys.argv[3]
  stokes_parm="I"
  for i in range(0,int(sys.argv[4]),5):
     lwimagerMS(msName, imageName, stokes_parm, i)

  stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
  print("Stop at %s" % stoptime)

  process_end = time.time()
  duration = (process_end - process_start)/3600.0
  print("Total run time: %7.2f hours" % duration)
