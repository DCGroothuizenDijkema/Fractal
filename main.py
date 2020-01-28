
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

def produce_mandelbrot_visualisation(example='one',show_fig=False,save_fig=True,file_dir='./out/',file_name='out.pdf'):

  try:
    os.makedirs(file_dir) 
  except OSError as err:
    if err.errno != errno.EEXIST:
      raise

  if example=='one':
    iterations=sample_mandelbrot((-0.7435669,-0.1314023),0.0022878,0.0022878,5001,5001,1000)
    plot_mandelbrot(iterations,file_name=file_dir+file_name)
    return

if __name__=='__main__':
  produce_mandelbrot_visualisation()
