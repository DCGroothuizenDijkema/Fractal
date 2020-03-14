
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
from matplotlib.colors import ListedColormap

from huygens.interf import c_matrix,c_vector,c_pointer

__all__=['sample_mandelbrot','plot_mandelbrot','sample_newton','plot_newton','plot_newton_roots','plot_newton_iteration']

# load the lib
_libc=ct.cdll.LoadLibrary('./bin/fractal.dll')

# extract the functions
_sample_mandelbrot=getattr(_libc,'?sample_mandelbrot@@YAHPEAPEAHHHHHNNNN_N@Z')
_sample_newton=getattr(_libc,'?sample_newton@@YAHPEAPEAN0PEAPEAHPEANHHHHHNNNN_N@Z')
_assign_roots=getattr(_libc,'?assign_roots@@YAXQEBQEAHQEBQEBN1QEBN2HHH@Z')

# assign arg and return types
_sample_mandelbrot.argtypes=[ct.POINTER(ct.POINTER(ct.c_int)),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_double
  ,ct.c_double,ct.c_double,ct.c_double,ct.c_bool]
_sample_mandelbrot.restype=ct.c_int
_sample_newton.argtypes=[ct.POINTER(ct.POINTER(ct.c_double)),ct.POINTER(ct.POINTER(ct.c_double)),ct.POINTER(ct.POINTER(ct.c_int))
  ,ct.POINTER(ct.c_double),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_double,ct.c_double,ct.c_double
  ,ct.c_double,ct.c_bool]
_sample_newton.restype=ct.c_int
_assign_roots.argtypes=[ct.POINTER(ct.POINTER(ct.c_int)),ct.POINTER(ct.POINTER(ct.c_double)),ct.POINTER(ct.POINTER(ct.c_double))
  ,ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.c_int,ct.c_int,ct.c_int]
_assign_roots.restype=None

def plot_mandelbrot(iterations,limit,log=True,show_fig=False,save_fig=True,file_name='mandelbrot.pdf',fig_inches=(12,12),dpi=1200
  ,color_map=None):

  _,ax=_plot_setup(fig_inches)

  if log:
    iterations=np.log(iterations)
    limit=np.log(limit)
  
  # black out where the limit could not be found (in the mandelbrot set)
  # and set color map
  masked_iterations=np.ma.masked_where(iterations==limit,iterations)
  if color_map is None:
    color_map=cm.Spectral_r
  color_map.set_bad(color='black')
  # produce the figure
  ax.imshow(masked_iterations,cmap=color_map)

  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def sample_mandelbrot(central_point,x_span,y_span,x_resolution,y_resolution,max_itr,num_threads=1,verbose=False):

  # input setup
  startx=central_point[0]-x_span/2.
  starty=central_point[1]-y_span/2.
  endx=central_point[0]+x_span/2.
  endy=central_point[1]+y_span/2.
  # array to store the iteration count for each pixel
  tmp,act=c_matrix(ct.c_int,y_resolution,x_resolution)

  # call the library function
  limit=_sample_mandelbrot(
    tmp,ct.c_int(max_itr),ct.c_int(num_threads),ct.c_int(x_resolution),ct.c_int(y_resolution)
    ,ct.c_double(startx),ct.c_double(endx),ct.c_double(starty),ct.c_double(endy),ct.c_bool(verbose)
  )

  del tmp
  # flip the rows because [startx,stary] is stored in [0,0]
  return np.flipud(np.ctypeslib.as_array(act)),limit

def plot_newton_roots(roots,show_fig=False,save_fig=True,file_name='newtons_fractal_roots.pdf',fig_inches=(12,12),dpi=1200
  ,color_map=None):
  
  _,ax=_plot_setup(fig_inches)

  # black out where no root could be found and set color map
  masked_roots=np.ma.masked_where(roots==-1,roots)
  if color_map is None:
    color_map=cm.viridis
  color_map.set_bad(color='black')
  # produce the figure
  ax.imshow(masked_roots,cmap=color_map)

  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def plot_newton_iteration(iterations,limit,log=True,show_fig=False,save_fig=True,file_name='newtons_fractal_iterations.pdf'
  ,fig_inches=(12,12),dpi=1200,color_map=None):

  _,ax=_plot_setup(fig_inches)
  
  if log:
    iterations=np.log(iterations)
    limit=np.log(limit)

  # black out where no root could be found and set color map
  masked_iterations=np.ma.masked_where(iterations==limit,iterations)
  if color_map is None:
    color_map=cm.Spectral_r
  color_map.set_bad(color='black')
  # produce the figure
  ax.imshow(masked_iterations,cmap=color_map)

  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def plot_newton(roots,iterations,limit,colors,log=True,show_fig=False,save_fig=True,file_name='mandelbrot.pdf',fig_inches=(12,12),dpi=1200):

  _,ax=_plot_setup(fig_inches)

  # find all unique roots by the index values, ignoring where no root was found
  unique_roots=np.unique(roots)
  unique_roots=np.delete(unique_roots,np.where(unique_roots==-1))

  if log:
    iterations=np.log(iterations)
    limit=np.log(limit)

  # mask where no roots could be found
  limit_bool=iterations!=limit
  # plot where no roots could be found as black
  no_root=np.ma.masked_array(iterations,limit_bool)
  ax.imshow(no_root,cmap=ListedColormap([0,0,0]))
  # for each root, mask it and plot it
  for itr,root in enumerate(unique_roots):
    masked_roots=np.ma.masked_array(iterations,roots!=root)
    ax.imshow(iterations*masked_roots,cmap=colors[itr])

  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def sample_newton(coeffs,central_point,x_span,y_span,x_resolution,y_resolution,max_itr,num_threads=1,verbose=False):

  # input setup
  startx=central_point[0]-x_span/2.
  starty=central_point[1]-y_span/2.
  endx=central_point[0]+x_span/2.
  endy=central_point[1]+y_span/2.
  
  # arrays to store the real part of the root approached and the imaginary part of the root approached
  tmp_re,act_re=c_matrix(ct.c_double,y_resolution,x_resolution)
  tmp_im,act_im=c_matrix(ct.c_double,y_resolution,x_resolution)
  # array to store the iteration count to get to the root
  tmp_itr,act_itr=c_matrix(ct.c_int,y_resolution,x_resolution)
  # coefficients of the polynomial in library form
  poly_coeffs=c_vector(ct.c_double,len(coeffs),coeffs)
  poly_degree=len(coeffs)-1

  # call the library function
  limit=_sample_newton(
    tmp_re,tmp_im,tmp_itr,poly_coeffs,ct.c_int(max_itr),ct.c_int(num_threads),ct.c_int(poly_degree),ct.c_int(x_resolution)
    ,ct.c_int(y_resolution),ct.c_double(startx),ct.c_double(endx),ct.c_double(starty),ct.c_double(endy),ct.c_bool(verbose)
  )

  del tmp_itr
  # convert roots to a numpy array
  roots=1j*np.ctypeslib.as_array(act_im)
  roots+=np.ctypeslib.as_array(act_re)

  # determine the actual roots
  actuals=np.roots(np.flip(coeffs))
  # real and imaginary parts of the actual roots
  roots_re=c_vector(ct.c_double,len(actuals),np.real(actuals))
  roots_im=c_vector(ct.c_double,len(actuals),np.imag(actuals))
  # array to store which root the approximation is
  tmp_ind,act_ind=c_matrix(ct.c_int,y_resolution,x_resolution)
  
  # call the library function
  _assign_roots(tmp_ind,tmp_re,tmp_im,roots_re,roots_im,ct.c_int(poly_degree),ct.c_int(x_resolution),ct.c_int(y_resolution))

  del tmp_re
  del tmp_im

  # flip the rows because [startx,stary] is stored in [0,0]
  return np.flipud(roots),np.flipud(np.ctypeslib.as_array(act_ind)),np.flipud(np.ctypeslib.as_array(act_itr)),limit

def _plot_setup(fig_inches):
  
  fig,ax=plt.subplots()
  # remove all borders
  fig.subplots_adjust(0,0,1,1)
  # remove all ticks and labels
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  # scale
  fig.set_size_inches(fig_inches)

  return fig,ax