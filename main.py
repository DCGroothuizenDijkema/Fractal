
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
#                                                                                                                                       //
# main.py                                                                                                                               //
#                                                                                                                                       //
# D. C. Groothuizen Dijkema - January, 2020                                                                                             //
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

# Main file to produce visualisations of fractals


import errno
import os

from fractal import *

def produce_mandelbrot_visualisation(example='zoom_level_zero',resolution=10001,show_fig=False,save_fig=True,file_dir='./out/'
  ,file_name='out.pdf'):

  try:
    os.makedirs(file_dir) 
  except OSError as err:
    if err.errno != errno.EEXIST:
      raise

  dx=3
  dy=2

  if example=='zoom_level_zero':
    span=1
    iterations,limit=sample_mandelbrot((-0.7,0.0),dx*span,dy*span,dx*resolution,dy*resolution,1000)
    plot_mandelbrot(iterations,limit,file_name=file_dir+file_name,fig_inches=(6*dx,6*dy))
  elif example=='zoom_level_one':
    span=0.007
    iterations,limit=sample_mandelbrot((-0.744030,0.1263),dx*span,dy*span,dx*resolution,dy*resolution,1000)
    plot_mandelbrot(iterations,limit,file_name=file_dir+file_name,fig_inches=(6*dx,6*dy),dpi=300)
  elif example=='zoom_level_two':
    span=0.0065
    iterations,limit=sample_mandelbrot((-0.775,0.121),dx*span,dy*span,dx*resolution,dy*resolution,1000)
    plot_mandelbrot(iterations,limit,file_name=file_dir+file_name,fig_inches=(6*dx,6*dy),dpi=300)
  elif example=='zoom_level_three':
    span=0.0011439
    iterations,limit=sample_mandelbrot((-0.7439668,0.1314023),dx*span,dy*span,dx*resolution,dy*resolution,1000)
    plot_mandelbrot(iterations,limit,file_name=file_dir+file_name,fig_inches=(6*dx,6*dy),dpi=1200)

if __name__=='__main__':
  produce_mandelbrot_visualisation('zoom_level_three')
