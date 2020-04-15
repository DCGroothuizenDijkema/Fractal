
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
#                                                                                                                                       //
# main.py                                                                                                                               //
#                                                                                                                                       //
# D. C. Groothuizen Dijkema - January, 2020                                                                                             //
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

# Main file to produce visualisations of fractals


from matplotlib.colors import LinearSegmentedColormap

from fractal import sample_mandelbrot,plot_mandelbrot,sample_newton,plot_newton

def produce_mandelbrot_visualisation(example='zoom_level_zero',fractal_resolution=3001,limit=1000,show_fig=False,save_fig=True
  ,file_name='mandelbrot.pdf',dpi=1200,num_threads=1,verbose=False):
  '''
  Calculate a portion of the Mandelbrot Set and produce a visualisation of it.
  
  Parameters
  ----------
  example : string, optional
    - The named example to produce
      One of: 'zoom_level_zero','zoom_level_one','zoom_level_two','zoom_level_three'
  fractal_resolution : int, optional
    - The number of pixels to use in the calculation
  limit : int, optional
    - The maximum number of iterations allowed
  show_fig : bool, optional
    - If the visualisation should be shown.
  save_fig : bool, optional
    - If the visualisation should be saved.
  file_name : string, optional
    - The name of the output.
  dpi : int, optional
    - Plot resolution.
  num_threads : int, optional
    - The number of threads to execute on.
  verbose : bool, optional.
    - For verbose output.

  '''
  # constant figure parameters
  dx=3
  dy=2

  # get example parameters
  if example=='zoom_level_zero':
    span=1
    centre=(-0.7,0.0)
  elif example=='zoom_level_one':
    span=0.007
    centre=(-0.744030,0.1263)
  elif example=='zoom_level_two':
    span=0.0065
    centre=(-0.775,0.121)
  elif example=='zoom_level_three':
    span=0.0011439
    centre=(-0.7439668,0.1314023)
  else:
    # invalid example
    examples=['zoom_level_zero','zoom_level_one','zoom_level_two','zoom_level_three']
    raise ValueError('`example` must be one of {}'.format(examples))
  
  # determine the fractal
  iterations,limit=sample_mandelbrot(centre,dx*span,dy*span,dx*fractal_resolution,dy*fractal_resolution,limit,num_threads=num_threads
    ,verbose=verbose)
  # produce the visualisation
  plot_mandelbrot(iterations,limit,file_name=file_name,fig_inches=(6*dx,6*dy),dpi=dpi,show_fig=show_fig,save_fig=save_fig)

def produce_newton_visualisation(example='cubic_zero',fractal_resolution=3001,limit=1000,show_fig=False,save_fig=True
  ,file_name='newton.pdf',dpi=1200,num_threads=1,verbose=False):
  '''
  Calculate a Newton's fractal and produce a visualisation of it.
  
  Parameters
  ----------
  example : string, optional
    - The named example to produce
      One of: 'cubic_zero','cubic_one','cubic_two','quartic_zero','quartic_one','quartic_two','pentic_zero','sextic_zero_zoom_level_zero',
        'sextic_zero_zoom_level_one','sextic_zero_zoom_level_two'
  fractal_resolution : int, optional
    - The number of pixels to use in the calculation
  limit : int, optional
    - The maximum number of iterations allowed
  show_fig : bool, optional
    - If the visualisation should be shown.
  save_fig : bool, optional
    - If the visualisation should be saved.
  file_name : string, optional
    - The name of the output.
  dpi : int, optional
    - Plot resolution.
  num_threads : int, optional
    - The number of threads to execute on.
  verbose : bool, optional.
    - For verbose output.

  '''
  # some linear colour scales for plotting
  colors=[
      LinearSegmentedColormap.from_list('custom_colormap',['#f3c8ea','#6f185d'])
      ,LinearSegmentedColormap.from_list('custom_colormap',['#eaf3c8','#28300a'])
      ,LinearSegmentedColormap.from_list('custom_colormap',['#b3c5ef','#18326f'])
      ,LinearSegmentedColormap.from_list('custom_colormap',['#f3c85d','#6f1833'])
      ,LinearSegmentedColormap.from_list('custom_colormap',['#f3d1c8','#6f2a18'])
      ,LinearSegmentedColormap.from_list('custom_colormap',['#c8eaf3','#185d6f'])
    ]

  # constant figure parameters
  centre=(0,0)
  span=2
  dx=2
  dy=2

  # get example parameters
  if example=='cubic_zero':
    coeffs=[-1,0,0,1]
  elif example=='cubic_one':
    coeffs=[2,-2,0,1]
  elif example=='cubic_two':
    coeffs=[4,-2,0,1]
  elif example=='quartic_zero':
    coeffs=[0,-1,0,0,3]
  elif example=='quartic_one':
    coeffs=[0,-2,2,-2,1]
  elif example=='quartic_two':
    coeffs=[2,-2,1,-2,2]
  elif example=='pentic_zero':
    coeffs=[-1,0,0,0,0,1]
  elif example=='sextic_zero_zoom_level_zero':
    coeffs=[1,2,8,3,5,0,1]
    # in this case, the geometry is more interesting with a rectangular view
    dx=4
    dy=3
  elif example=='sextic_zero_zoom_level_one':
    coeffs=[1,2,8,3,5,0,1]
    # this is a zoom
    centre=(0.45,-1.4)
    span=0.3
  elif example=='sextic_zero_zoom_level_two':
    coeffs=[1,2,8,3,5,0,1]
    # this is a zoom
    centre=(0.535,-1.24)
    span=0.08
  else:
    # invalid example
    examples=['cubic_zero','cubic_one','cubic_two','quartic_zero','quartic_one','quartic_two','pentic_zero','sextic_zero_zoom_level_zero'
      ,'sextic_zero_zoom_level_one','sextic_zero_zoom_level_two']
    raise ValueError('`example` must be one of {}'.format(examples))

  # determine the fractal
  _,idx,itr,limit=sample_newton(coeffs,centre,dx*span,dy*span,dx*fractal_resolution,dy*fractal_resolution,limit,num_threads,verbose)
  # produce the visualisation
  plot_newton(idx,itr,limit,colors,file_name=file_name,fig_inches=(6*dx,6*dy),dpi=dpi,show_fig=show_fig,save_fig=save_fig)

if __name__=='__main__':
  produce_mandelbrot_visualisation()
  produce_newton_visualisation()
