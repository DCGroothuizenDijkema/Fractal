
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//
#                                                                                                                                       //
# main.py                                                                                                                               //
#                                                                                                                                       //
# D. C. Groothuizen Dijkema - January, 2020                                                                                             //
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+//

# Main file to produce visualisations of fractals


from fractal import sample_mandelbrot,plot_mandelbrot

def produce_mandelbrot_visualisation(example='zoom_level_zero',fractal_resolution=3001,show_fig=False,save_fig=True,
  file_name='mandelbrot.pdf',dpi=1200,verbose=False):

  dx=3
  dy=2

  if example=='zoom_level_zero':
    span=1
    iterations,limit=sample_mandelbrot((-0.7,0.0),dx*span,dy*span,dx*fractal_resolution,dy*fractal_resolution,1000,verbose=verbose)
    plot_mandelbrot(iterations,limit,file_name=file_name,fig_inches=(6*dx,6*dy),dpi=dpi)
  elif example=='zoom_level_one':
    span=0.007
    iterations,limit=sample_mandelbrot((-0.744030,0.1263),dx*span,dy*span,dx*fractal_resolution,dy*fractal_resolution,1000,verbose=verbose)
    plot_mandelbrot(iterations,limit,file_name=file_name,fig_inches=(6*dx,6*dy),dpi=dpi)
  elif example=='zoom_level_two':
    span=0.0065
    iterations,limit=sample_mandelbrot((-0.775,0.121),dx*span,dy*span,dx*fractal_resolution,dy*fractal_resolution,1000,verbose=verbose)
    plot_mandelbrot(iterations,limit,file_name=file_name,fig_inches=(6*dx,6*dy),dpi=dpi)
  elif example=='zoom_level_three':
    span=0.0011439
    iterations,limit=sample_mandelbrot((-0.7439668,0.1314023),dx*span,dy*span,dx*fractal_resolution,dy*fractal_resolution,1000
      ,verbose=verbose)
    plot_mandelbrot(iterations,limit,file_name=file_name,fig_inches=(6*dx,6*dy),dpi=dpi)

if __name__=='__main__':
  produce_mandelbrot_visualisation()
