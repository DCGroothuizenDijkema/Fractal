
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

def produce_mandelbrot_visualisation(example='zoom_level_zero',resolution=1001,show_fig=False,save_fig=True,file_dir='./out/'
  ,file_name='out.pdf'):

  try:
    os.makedirs(file_dir) 
  except OSError as err:
    if err.errno != errno.EEXIST:
      raise

  if example=='zoom_level_zero':
    dx=3
    dy=2
    iterations,limit=sample_mandelbrot((-0.7,0.0),3,2,dx*resolution,dy*resolution,1000)
    plot_mandelbrot(iterations,limit,file_name=file_dir+file_name,fig_inches=(6*dx,6*dy))

if __name__=='__main__':
  produce_mandelbrot_visualisation('zoom_level_zero')
