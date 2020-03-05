
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
from matplotlib.colors import LinearSegmentedColormap,ListedColormap

from huygens.interf import c_matrix,c_vector,c_pointer

__all__=['sample_mandelbrot','plot_mandelbrot','sample_newton','plot_newton','plot_newton_roots','plot_newton_iteration']

# load the lib
_libc=ct.cdll.LoadLibrary('./bin/fractal.dll')

# extract the functions
_sample_mandelbrot=getattr(_libc,'?sample_mandelbrot@@YAXPEAPEAHHHHHQEAHNNNN_N@Z')
_sample_newton=getattr(_libc,'?sample_newton@@YAXPEAPEAN0PEAPEAHPEANHHHHHQEAHNNNN_N@Z')
_assign_roots=getattr(_libc,'?assign_roots@@YAXQEBQEAHQEBQEBN1QEBN2HHH@Z')

# assign arg and return types
_sample_mandelbrot.argtypes=[ct.POINTER(ct.POINTER(ct.c_int)),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.POINTER(ct.c_int),ct.c_double
  ,ct.c_double,ct.c_double,ct.c_double,ct.c_bool]
_sample_mandelbrot.restype=None
_sample_newton.argtypes=[ct.POINTER(ct.POINTER(ct.c_double)),ct.POINTER(ct.POINTER(ct.c_double)),ct.POINTER(ct.POINTER(ct.c_int))
  ,ct.POINTER(ct.c_double),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.POINTER(ct.c_int),ct.c_double,ct.c_double,ct.c_double
  ,ct.c_double,ct.c_bool]
_sample_newton.restype=None
_assign_roots.argtypes=[ct.POINTER(ct.POINTER(ct.c_int)),ct.POINTER(ct.POINTER(ct.c_double)),ct.POINTER(ct.POINTER(ct.c_double))
  ,ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.c_int,ct.c_int,ct.c_int]
_assign_roots.restype=None

def plot_mandelbrot(iterations,limit,log=True,show_fig=False,save_fig=True,file_name='mandelbrot.pdf',fig_inches=(12,12),dpi=1200
  ,color_map=None):
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
  if color_map is None:
    color_map=cm.Spectral_r
  color_map.set_bad(color='black')

  ax.imshow(masked_iterations,cmap=color_map)
  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def sample_mandelbrot(central_point,x_span,y_span,x_resolution,y_resolution,max_itr,num_threads=1,verbose=False):
  startx=central_point[0]-x_span/2.
  starty=central_point[1]-y_span/2.
  endx=central_point[0]+x_span/2.
  endy=central_point[1]+y_span/2.
  limit=c_pointer(ct.c_int,0)
  
  tmp,act=c_matrix(ct.c_int,y_resolution,x_resolution)
  _sample_mandelbrot(
    tmp,ct.c_int(max_itr),ct.c_int(num_threads),ct.c_int(x_resolution),ct.c_int(y_resolution),limit
    ,ct.c_double(startx),ct.c_double(endx),ct.c_double(starty),ct.c_double(endy),ct.c_bool(verbose)
  )

  del tmp
  return np.flipud(np.ctypeslib.as_array(act)),limit.contents.value

def plot_newton_roots(roots,show_fig=False,save_fig=True,file_name='newtons_fractal_roots.pdf',fig_inches=(12,12),dpi=1200
  ,color_map=None):
  fig,ax=plt.subplots()
  fig.subplots_adjust(0,0,1,1)

  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  
  fig.set_size_inches(fig_inches)

  masked_roots=np.ma.masked_where(roots==-1,roots)
  if color_map is None:
    color_map=cm.viridis
  color_map.set_bad(color='black')

  ax.imshow(masked_roots,cmap=color_map)
  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def plot_newton_iteration(iterations,limit,log=True,show_fig=False,save_fig=True,file_name='newtons_fractal_iterations.pdf'
  ,fig_inches=(12,12),dpi=1200,color_map=None):
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
  if color_map is None:
    color_map=cm.Spectral_r
  color_map.set_bad(color='black')

  ax.imshow(masked_iterations,cmap=color_map)
  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def plot_newton(roots,iterations,limit,colors,log=True,show_fig=False,save_fig=True,file_name='mandelbrot.pdf',fig_inches=(12,12),dpi=1200):
  fig,ax=plt.subplots()
  fig.subplots_adjust(0,0,1,1)

  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  
  fig.set_size_inches(fig_inches)

  unique_roots=np.unique(roots)
  unique_roots=np.delete(unique_roots,np.where(unique_roots==-1))

  limit_index=np.where(iterations==limit)
  limit_bool=iterations!=limit
  iterations=np.log(iterations)
  
  v1b=np.ma.masked_array(iterations,roots!=0)
  v1c=np.ma.masked_array(iterations,roots!=1)
  v1d=np.ma.masked_array(iterations,roots!=2)
  v1a=np.ma.masked_array(iterations,limit_bool)

  ax.imshow(v1a,cmap=ListedColormap([0,0,0]))
  ax.imshow(iterations*v1b,cmap=LinearSegmentedColormap.from_list('cust1',['#f3c8ea','#6f185d']))
  ax.imshow(iterations*v1c,cmap=LinearSegmentedColormap.from_list('cust2',['#eaf3c8','#5d6f18']))
  ax.imshow(iterations*v1d,cmap=LinearSegmentedColormap.from_list('cust3',['#c8eaf3','#06161b']))

  if show_fig:
    plt.show()
  if save_fig:
    plt.savefig(file_name,dpi=dpi)

def sample_newton(coeffs,central_point,x_span,y_span,x_resolution,y_resolution,max_itr,num_threads=1,verbose=False):
  startx=central_point[0]-x_span/2.
  starty=central_point[1]-y_span/2.
  endx=central_point[0]+x_span/2.
  endy=central_point[1]+y_span/2.
  limit=c_pointer(ct.c_int,0)
  
  tmp_re,act_re=c_matrix(ct.c_double,y_resolution,x_resolution)
  tmp_im,act_im=c_matrix(ct.c_double,y_resolution,x_resolution)
  tmp_itr,act_itr=c_matrix(ct.c_int,y_resolution,x_resolution)

  poly_coeffs=c_vector(ct.c_double,len(coeffs),coeffs)
  poly_degree=len(coeffs)-1

  _sample_newton(
    tmp_re,tmp_im,tmp_itr,poly_coeffs,ct.c_int(max_itr),ct.c_int(num_threads),ct.c_int(poly_degree)
    ,ct.c_int(x_resolution),ct.c_int(y_resolution),limit,ct.c_double(startx),ct.c_double(endx)
    ,ct.c_double(starty),ct.c_double(endy),ct.c_bool(verbose)
  )

  del tmp_itr

  roots=1j*np.ctypeslib.as_array(act_im)
  roots+=np.ctypeslib.as_array(act_re)

  actuals=np.roots(coeffs)
  tmp_ind,act_ind=c_matrix(ct.c_int,y_resolution,x_resolution)
  roots_re=c_vector(ct.c_double,len(actuals),np.real(actuals))
  roots_im=c_vector(ct.c_double,len(actuals),np.imag(actuals))

  _assign_roots(tmp_ind,tmp_re,tmp_im,roots_re,roots_im,ct.c_int(poly_degree),ct.c_int(x_resolution),ct.c_int(y_resolution))

  del tmp_re
  del tmp_im

  return np.flipud(roots),np.flipud(np.ctypeslib.as_array(act_ind)),np.flipud(np.ctypeslib.as_array(act_itr)),limit.contents.value
