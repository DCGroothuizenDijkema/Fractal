
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
import matplotlib.cm as cm

from huygens.interf import c_matrix,c_pointer

__all__=['sample_mandelbrot','plot_mandelbrot']

# load the lib
_libc=ct.cdll.LoadLibrary('./bin/fractal.dll')

# extract the functions
_sample_mandelbrot=getattr(_libc,'?sample_mandelbrot@@YAXPEAPEAHHHHQEAHNNNN@Z')

# assign arg and return types
_sample_mandelbrot.argtypes=[ct.POINTER(ct.POINTER(ct.c_int)),ct.c_int,ct.c_int,ct.c_int,ct.POINTER(ct.c_int),ct.c_double,ct.c_double,ct.c_double,ct.c_double]
_sample_mandelbrot.restype=None

def plot_mandelbrot(iterations,limit,log=True,file_name='out.png',fig_inches=(12,12),dpi=1200):
  fig,ax=plt.subplots()
  fig.subplots_adjust(0,0,1,1)

  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xticklabels([])
  ax.set_yticklabels([])

  if log:
    iterations=np.log(iterations)
    limit=np.log(limit)
  
  fig.set_size_inches(fig_inches)

  masked_iterations=np.ma.masked_where(iterations==limit,iterations)
  cmap=cm.Spectral_r
  cmap.set_bad(color='black')

  ax.imshow(masked_iterations,cmap=cmap)
  plt.savefig(file_name,dpi=dpi)

def sample_mandelbrot(central_point,x_span,y_span,x_resolution,y_resolution,max_itr):
  startx=central_point[0]-x_span/2.
  starty=central_point[1]-y_span/2.
  endx=central_point[0]+x_span/2.
  endy=central_point[1]+y_span/2.
  limit=c_pointer(ct.c_int,0)
  
  tmp,act=c_matrix(ct.c_int,y_resolution,x_resolution)
  _sample_mandelbrot(
    tmp,ct.c_int(max_itr),ct.c_int(x_resolution),ct.c_int(y_resolution),limit
    ,ct.c_double(startx),ct.c_double(endx),ct.c_double(starty),ct.c_double(endy)
  )

  del tmp
  return np.ctypeslib.as_array(act),limit.contents.value
