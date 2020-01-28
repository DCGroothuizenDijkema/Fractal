
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
#                                                                                                                                       //
# fractal.py                                                                                                                            //
#                                                                                                                                       //
# D. C. Groothuizen Dijkema - January, 2020                                                                                             //
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

# Helper function file to produce visualisations of fractals


import ctypes as ct
import numpy as np
import matplotlib.pyplot as plt

from huygens.interf import c_matrix

__all__=['sample_mandelbrot','plot_mandelbrot']

# load the lib
_libc=ct.cdll.LoadLibrary('./bin/fractal.dll')

# extract the functions
_sample_mandelbrot=getattr(_libc,'?sample_mandelbrot@@YAXPEAPEAHHHHNNNN@Z')

# assign arg and return types
_sample_mandelbrot.argtypes=[ct.POINTER(ct.POINTER(ct.c_int)),ct.c_int,ct.c_int,ct.c_int,ct.c_double,ct.c_double,ct.c_double,ct.c_double]
_sample_mandelbrot.restype=None

def plot_mandelbrot(iterations,log=True,file_name='out.png',dpi=1200):
  fig,ax=plt.subplots()
  fig.subplots_adjust(0,0,1,1)

  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  if log:
    iterations=np.log(iterations)
  ax.imshow(iterations,cmap='Spectral_r')
  plt.savefig(file_name,dpi=dpi)

def sample_mandelbrot(central_point,x_span,y_span,x_resolution,y_resolution,max_itr):
  startx=central_point[0]-x_span/2.
  starty=central_point[1]-y_span/2.
  endx=central_point[0]+x_span/2.
  endy=central_point[1]+y_span/2.
  
  tmp,act=c_matrix(ct.c_int,y_resolution,x_resolution)
  _sample_mandelbrot(
    tmp,ct.c_int(max_itr),ct.c_int(x_resolution),ct.c_int(y_resolution)
    ,ct.c_double(startx),ct.c_double(endx),ct.c_double(starty),ct.c_double(endy)
  )

  del tmp
  return np.ctypeslib.as_array(act)
