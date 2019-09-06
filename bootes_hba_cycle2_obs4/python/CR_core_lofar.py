import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pl
from plot_util import *

import scipy
import pyfits as pf
import numpy as np
import scipy.signal as sps
from scipy import ndimage
import os
import time
import sys
try:
    from lofar import bdsm
except ImportError:
    print "no lofar"
import pywcs as pw
#6from spam import *
#from wwpy.multi_proc import run_in_separate_process
#from wwpy.extraction import *
from fits_util import *
#import pylab as pl
from scipy.stats import pareto
import warnings
warnings.filterwarnings('ignore')
#from multi import run_in_separate_process

from matplotlib.ticker import MultipleLocator, LogFormatter, LogLocator, FuncFormatter
#mpl.rc_file('/home/wwilliams/.config/matplotlib/matplotlibrc')  # <-- the file containing your settings
#mpl.rc_file('/net/laak/data2/wwilliams/phd/gmrt_bootes/a2256_mosaic/fits/matplotlibrc') 
#mpl.rc_file('/home/wwilliams/.config/matplotlib/matplotlibrc')  # <-- the file containing your settings

from sky_util import angular_separation



# get source finding routines 
try:
    from find_sources_pybdsm import *
except ImportError:
    print " no find_sources"


def fig_save_many(f, name, types=[".png",".pdf"], dpi=200):
    for ext in types:
        f.savefig(name+ext, dpi=dpi)
    return

#def paper_single(TW = 6.64, AR = 0.74, FF = 1.):
    ##import matplotlib as mpl
    ## textwidht = 42pc = 42 * 12 pt = 42 * 12 * 1/72.27 inches
    ## columnsep = 2pc
    ## ... colwidth = 20pc
    #'''
    #TW = 3.32
    #AR = 0.74
    #FF = 1.
    ##mpl.rc('figure', figsize=(4.5,3.34), dpi=200)
    #mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    #mpl.rc('figure.subplot', left=0.18, right=0.97, bottom=0.18, top=0.9)
    #mpl.rc('lines', linewidth=1.0, markersize=4.0)
    #mpl.rc('font', size=9.0, family="serif", serif="CM")
    #mpl.rc('xtick', labelsize='small')
    #mpl.rc('ytick', labelsize='small')
    #mpl.rc('axes', linewidth=0.75)
    #mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    #mpl.rc('text', usetex=True) 
    #mpl.rc('savefig', dpi=300)
    #'''
    
    ##mpl.rc('figure', figsize=(4.5,3.34), dpi=200)
    #mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=100)
    #mpl.rc('figure.subplot', left=0.15, right=0.95, bottom=0.15, top=0.92)
    #mpl.rc('lines', linewidth=1.75, markersize=8.0, markeredgewidth=0.75)
    #mpl.rc('font', size=18.0, family="serif", serif="CM")
    ##mpl.rc('font', size=18.0, serif="CM")
    #mpl.rc('xtick', labelsize='small')
    #mpl.rc('ytick', labelsize='small')
    #mpl.rc('xtick.major', width=1.0, size=8)
    #mpl.rc('ytick.major', width=1.0, size=8)
    #mpl.rc('xtick.minor', width=1.0, size=4)
    #mpl.rc('ytick.minor', width=1.0, size=4)
    #mpl.rc('axes', linewidth=1.5)
    #mpl.rc('legend', fontsize='small', numpoints=1, labelspacing=0.4, frameon=False) 
    #mpl.rc('text', usetex=True) 
    #mpl.rc('savefig', dpi=300)
    
    
    
#def paper_single_mult_ax(nrows=1, ncols=1, **kwargs):
    ##import matplotlib as mpl
    #paper_single(FF=max(nrows,ncols))
    #f, ax = pl.subplots(nrows=nrows, ncols=ncols, **kwargs)
    #pl.minorticks_on()
    #ylocator6 = pl.MaxNLocator(5)
    #xlocator6 = pl.MaxNLocator(6)
    #if len(ax.shape) > 1:
        #for axrow in ax:
            #for axcol in axrow:
                #axcol.xaxis.set_major_locator(xlocator6)
                #axcol.yaxis.set_major_locator(ylocator6)
    #else:
        #for axcol in ax:
            #axcol.xaxis.set_major_locator(xlocator6)
            #axcol.yaxis.set_major_locator(ylocator6)
    #return f, ax
    
    
#def paper_single_ax():
    ##import matplotlib as mpl
    #paper_single()
    #f = pl.figure()
    #ax = pl.subplot(111)
    #pl.minorticks_on()
    #ylocator6 = pl.MaxNLocator(5)
    #xlocator6 = pl.MaxNLocator(6)
    #ax.xaxis.set_major_locator(xlocator6)
    #ax.yaxis.set_major_locator(ylocator6)
    #return f, ax


#def paper_double_ax():
    ##import matplotlib as mpl
    #paper_single(TW = 12)
    #f = pl.figure()
    #ax = pl.subplot(111)
    #pl.minorticks_on()
    #ylocator6 = pl.MaxNLocator(5)
    #xlocator6 = pl.MaxNLocator(6)
    #ax.xaxis.set_major_locator(xlocator6)
    #ax.yaxis.set_major_locator(ylocator6)
    #return f, ax

#def paper_double_mult_ax(nrows=1, ncols=1, setticks=True, **kwargs):
    ##import matplotlib as mpl
    #paper_single()
    ##TW = 6.97
    ##AR = 0.74
    ##FF = 1.
    ##mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    ##mpl.rc('figure.subplot', left=0.08, right=0.95, bottom=0.08, top=0.95)
    ##mpl.rc('font', size=12.0, family="serif", serif="CM")
    ##f = pl.figure()
    ##ax = pl.subplot(111)
    ##pl.minorticks_on()
    ##ylocator6 = pl.MaxNLocator(7)
    ##xlocator6 = pl.MaxNLocator(8)
    ##ax.xaxis.set_major_locator(xlocator6)
    ##ax.yaxis.set_major_locator(ylocator6)
    
    
    #TW = 6.97*2
    #AR = 0.74
    #FF = 1.
    #mpl.rc('figure', figsize=(FF*TW, FF*TW*AR), dpi=200)
    #mpl.rc('figure.subplot', left=0.1, right=0.97, bottom=0.1, top=0.97)
    #mpl.rc('font', size=24.0, family="serif", serif="CM")
    ##mpl.rc('font', size=24.0, serif="CM")
    
    #f, ax = pl.subplots(nrows=nrows, ncols=ncols, **kwargs)
    #pl.minorticks_on()
    #if setticks:
        #ylocator6 = pl.MaxNLocator(5)
        #xlocator6 = pl.MaxNLocator(6)
        #if len(ax.shape) > 1:
            #for axrow in ax:
                #for axcol in axrow:
                    #axcol.xaxis.set_major_locator(xlocator6)
                    #axcol.yaxis.set_major_locator(ylocator6)
        #else:
            #for axcol in ax:
                #axcol.xaxis.set_major_locator(xlocator6)
                #axcol.yaxis.set_major_locator(ylocator6)
    #return f, ax



def makeGaussian(size,amplitude, x0,y0,dx,dy,rota , fwhm = 3, center=None):
    """ Make a square gaussian kernel.
 
size is the length of a side of the square
fwhm is full-width-half-maximum, which
can be thought of as an effective radius.
"""
 
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
        
    rota = np.pi/180. * rota
    
    x0p = x0 * np.cos(rota) - y0 * np.sin(rota)
    y0p = x0 * np.sin(rota) + y0 * np.cos(rota)
    
    xp = x * np.cos(rota) - y * np.sin(rota)
    yp = x * np.sin(rota) + y * np.cos(rota)
        
    return amplitude*np.exp( -( (((x0p-xp)**2.)/(dx**2.)) + (((y0p-yp)**2.)/(dy**2.)) ) )

def twodgaussian(x,y,amplitude, x0,y0,dx,dy,rota):

    rota = np.pi/180. * rota
    x0p = x0 * np.cos(rota) - y0 * np.sin(rota)
    y0p = x0 * np.sin(rota) + y0 * np.cos(rota)
            
    xp = x * np.cos(rota) - y * np.sin(rota)
    yp = x * np.sin(rota) + y * np.cos(rota)
    g = amplitude*np.exp( -( (((x0p-xp)**2.)/(dx**2.)) + (((y0p-yp)**2.)/(dy**2.)) ) )
    
    return g


def add_gaussian(data,Speak,x0,y0,dx,dy,bpa):
  #print flux,x0,y0,dx,dy,bpa
  #Speak = flux/(2.*np.pi*dx*dy)
  #result = True
  #ymin = 0
  C = 6
  if dx/dy < 0.2:
    C = 8
  #xmax,ymax =  data.shape
  size = np.ceil(max(C*dx, C*dy))
  x1 = np.floor(x0) - size
  x2 = np.ceil(x0) + size
  y1 = np.floor(y0) - size
  y2 = np.ceil(y0) + size
  #print x1,x2
  #print y1,y2
  #for xi in range( x1, x2 ):
    #for yi in range( y1, y2 ):
      ##print xi,yi
      #data[xi,yi] += Speak*twodgaussian(xi,yi,1,x0,y0,dx,dy,bpa)
      #if rms[xi,y1] == 0.:
	#result = True
	#return data, result
    
  x = np.arange(-size, size, 1, float)
  y = x[:,np.newaxis]
  rota = np.pi/180. * bpa
  x0p = x0 * np.cos(rota) - y0 * np.sin(rota)
  y0p = x0 * np.sin(rota) + y0 * np.cos(rota)
            
  xp = x * np.cos(rota) - y * np.sin(rota)
  yp = x * np.sin(rota) + y * np.cos(rota)
  g = np.exp( -( (((x0p-xp)**2.)/(dx**2.)) + (((y0p-yp)**2.)/(dy**2.)) ) )
  
  #print size, x1,x2,y1,y2, data[x1:x2,y1:y2].shape, g.shape, Speak
  data[x1:x2,y1:y2] += Speak*g
    
  return data
  
def add_gaussian_point(data,Speak,x0,y0,dx,gaus):
  #print flux,x0,y0,dx,dy,bpa
  #Speak = flux/(2.*np.pi*dx*dy)
  #result = True
  #ymin = 0  
  #print size, x1,x2,y1,y2, data[x1:x2,y1:y2].shape, g.shape, Speak
  
  C = 6
  #xmax,ymax =  data.shape
  size = np.ceil(max(C*dx, C*dx))
  x1 = np.floor(x0) - size
  x2 = np.ceil(x0) + size
  y1 = np.floor(y0) - size
  y2 = np.ceil(y0) + size
  
  data[x1:x2,y1:y2] += Speak*gaus
    
  return data


def random_flux(Fmin=1e-3, Fmax=10.0, slope=-0.8):
  # fmin and fmax in Jy
  #F = pareto.rvs(0.8)*10.*sig
  F = pareto.rvs(-1*slope,loc=-Fmin, scale=Fmin, size=1) #*5.*sig
  satisfy = False
  #print F
  while not satisfy:
    #print 'repeat', F, F < 100.*sig
    if (F < Fmax) and (F > Fmin):
      satisfy = True
    else:
      F = pareto.rvs(-1*slope,loc=-Fmin, scale=Fmin, size=1)#*0.1*sig
  return F # + 0.1*sig


#def random_flux_powerlaw(alpha=-1.6, fmin=0.1, fmax=10.):
  #found = False
  #while not found:
    #x = np.random.triangular(fmin,fmin,fmax)
    ##x = np.random.uniform(fmin,fmax)
    #y = np.random.uniform(10**(alpha*np.log10(fmin)), 10**(alpha*np.log10(fmax))) 
    #if y < 10**(alpha*np.log10(x)):
      #frand = x
      #found = True
  #return frand

def random_flux_powerlaw(alpha=-0.6, Fmin=0.01, Fmax=10.):
  found = False
  while not found:
    x = np.random.uniform(Fmin,Fmax)
    y = np.random.uniform(0, Fmin**(alpha)) 
    if y < x**(alpha):
      frand = x
      found = True
  return frand
  
def gaussian_circ(base, height, width):
    """Returns a gaussian function with the given parameters"""
    return lambda r: base - height*np.exp(-((r/width)**2/2))
  
def psf(shape, sig_beam_maj, sig_beam_min, beam_bpa, C=10):
    #make a convolved_point to scale with peak
   
    colv_point  = np.zeros(shape)
    x0 = shape[0]/2
    y0 = shape[1]/2
    Theta = beam_bpa*np.pi/180.  # in rad
    #Theta = 0.
    for xi in range( 0, shape[0] ):
        for yi in range( 0, shape[1] ):
            #print xi,yi
            x2 = np.cos(Theta)*(xi-x0)-np.sin(Theta)*(yi-y0)#+x0
            y2 = np.sin(Theta)*(xi-x0)+np.cos(Theta)*(yi-y0)#+y0
            colv_point[xi,yi] += np.exp( -( (((x2)**2.)/(2.*sig_beam_maj**2.)) + (((y2)**2.)/(2.*sig_beam_min**2.)) ) )
            #colv_point[xi,yi] += np.exp( -( (((xi-x0)**2.)/(2.*sig_beam_maj**2.)) + (((yi-y0)**2.)/(2.*sig_beam_min**2.)) ) )
    return colv_point
    

def make_model_flatmap(model_map_name, pb_map, model_flatmap_name):
    
    head = pf.getheader(model_map_name)
    data = pf.getdata(model_map_name)
    pb = pf.getdata(pb_map)
    
    flatdata = data*pb
    
    pf.writeto(model_flatmap_name, flatdata, header=head, clobber=True)
    
    return


def make_model_map(model_map_name, model_cat, residmap, sourcemap, Fmin, Fmax, alpha, Nsrc = 1000, Fracpoint=1., verbose=3):

  print "making ",model_map_name
  
  if isinstance(Nsrc, np.float):
      Nsrc = int(Nsrc)
    
  os.system('cp %s %s.old' %(residmap,residmap))
  
  ## get some info from the source image ##
  sourcehead = pf.getheader(sourcemap)
  beam_maj = sourcehead.get('BMAJ')
  beam_min = sourcehead.get('BMIN')
  beam_pa = sourcehead.get('BPA')
    # or it is not there, check the history
  if beam_maj is None:
      hist = sourcehead.get_history()
      for i in range(len(hist)):
            t = hist[i]
            if 'BMAJ=' in t:
                tt = t.split()
                beam_maj = float(tt[tt.index('BMAJ=')+1])
            if 'BMIN=' in t:
                tt = t.split()
                beam_min = float(tt[tt.index('BMIN=')+1])
            if 'BPA=' in t:
                tt = t.split()
                beam_pa = float(tt[tt.index('BPA=')+1])
  mapwcs = pw.WCS(sourcehead)
  
  ps = abs(sourcehead.get('CDELT1'))*3600.  #arcsec per pixel
  fwhm_beam_pix = beam_maj*3600./ps    #pixels
  fwhm_beam_maj_pix = beam_maj*3600./ps    #pixels
  fwhm_beam_min_pix = beam_min*3600./ps    #pixels
  fwhm_beam = beam_maj*3600.    #arcsec
  fwhm_beam_maj = beam_maj*3600.    #arcsec
  fwhm_beam_min = beam_min*3600.    #arcsec
  
  ## get the residual data ##
  residall = pf.getdata(residmap)#.transpose()
  S = residall.shape
  Nax = len(S)
  if Nax > 2:
      residall = residall[0,0,:,:]
      
  resid=residall

  S = resid.shape
  Nax = len(S)
  model = np.zeros(resid.shape)

  
  modelfits = model_map_name
  #if os.path.isfile( modelfits ):
    #os.system(  'rm '+modelfits )
  #nx = S[0]
  #ny = S[1]
  #inhead['NAXIS1'] = int(nx)
  #inhead['NAXIS2'] = int(ny)
  #model1 = np.array([[model]])
  #hdu = pf.PrimaryHDU(model1, inhead)
  #hdu.writeto(modelfits)
  if verbose > 2: print 'shape:',model.shape
  
  
  fwhm_f = 2.*np.sqrt(2.*np.log(2.))
  sig_beam = fwhm_beam_pix/fwhm_f    ##pixels/
  sig_beam_maj = fwhm_beam_maj_pix/fwhm_f    ##pixels/
  sig_beam_min = fwhm_beam_min_pix/fwhm_f    ##pixels/
  
  Speak = np.zeros(Nsrc)
  Flux = np.zeros(Nsrc)
  bmaj = np.zeros(Nsrc)
  bmin = np.zeros(Nsrc)
  bpa = np.zeros(Nsrc)
  ra = np.zeros(Nsrc)
  dec = np.zeros(Nsrc)
  pnt = np.zeros(Nsrc)
  
  
  n_ext2 = 0
  if verbose > 1: print 'beam = %.3f pix' %(fwhm_beam_pix)
  if verbose > 1: print 'beam maj = %.3f pix' %(fwhm_beam_maj_pix)
  if verbose > 1: print 'beam min = %.3f pix' %(fwhm_beam_min_pix)
  if verbose > 1: print 'beam = %.3f arcsec' %(fwhm_beam)
  
  
  ##make a convolved_point to scale with peak
  #C = 6
  #int_beam = int(round(C*sig_beam,0))
  #colv_point  = np.zeros((2*int_beam+1,2*int_beam+1))
  #x0 = int_beam
  #y0 = int_beam
  #print 2*int_beam
  #for xi in range( 0, 2*int_beam+1 ):
    #for yi in range( 0, 2*int_beam+1 ):
      #colv_point[xi,yi] += twodgaussian(xi,yi,1., x0,y0,np.sqrt(2.)*sig_beam,np.sqrt(2.)*sig_beam,0.)
      
  #m_valx, m_valy =  colv_point.shape
  #print 'point kernel, %i x %i '  %(m_valx, m_valy)
  
  if verbose > 1: print 'sig_beam', sig_beam
  
  C = 6
  #xmax,ymax =  data.shape
  Psize = np.ceil(max(C*sig_beam, C*sig_beam))
  x0 = Psize
  y0 = Psize  
  xg = np.arange(0, 2*Psize, 1, float)
  yg = xg[:,np.newaxis]
  rota = np.pi/180. * (beam_pa-90)
  x0p = x0 * np.cos(rota) - y0 * np.sin(rota)
  y0p = x0 * np.sin(rota) + y0 * np.cos(rota)
  xp = xg * np.cos(rota) - yg * np.sin(rota)
  yp = xg * np.sin(rota) + yg * np.cos(rota)
  gaus_point = np.exp( -( (((x0p-xp)**2.)/(2.*sig_beam_maj**2.)) + (((y0p-yp)**2.)/(2.*sig_beam_min**2.)) ) )
  

 
  if verbose > 1: print 'gaussian point shape: ',gaus_point.shape
  
  #valid_coords    
  #x_good,y_good = np.where(np.isfinite(resid))
  nx,ny=resid.shape
  # set edges to nan
  resid[0,:] *= np.nan
  resid[-1,:] *= np.nan
  resid[:,0] *= np.nan
  resid[:,-1] *= np.nan
  badmask = np.isnan(resid)
  
  #import pylab as pl
  #pl.figure()
  #pl.imshow(badmask)
  badmask = ndimage.morphology.binary_dilation(badmask, iterations=int(1.5*Psize))  # grow mask in all directions by 50 pixels
  
  #pl.figure()
  #pl.imshow(badmask)
  #pl.show()
  x_good,y_good = np.where(~badmask)
  N_good = len(x_good)
  #i_cs = range(N_good)
  #np.random.shuffle(i_cs)
  
  xlist = 999.*np.ones(Nsrc)
  ylist = 999.*np.ones(Nsrc)
  # generate sources #
  for i in range(Nsrc):
    if i%100 == 0:
        print '%i of %i' %(i,Nsrc)
      
    #T1 = time.time()
      
    # decide if it is a point
    randpoint = np.random.uniform(0.,1.)
    if randpoint < Fracpoint:
        point = True
    else:
        point = False
        
    # repeat until we have a valid source (invalids occur on NaNs and near edges)
    result = False
    indi = 0
    nx , ny = model.shape
    #print  ' * %f mus' %((time.time()-T1)*1e6)
    while not result:
      result = True
      indi += 1
      N_good = len(x_good)
      i_c = np.random.randint(0,N_good)
      #i_c = i_cs[indi]
      #x = np.random.uniform(0,nx)
      #y = np.random.uniform(0,ny)
      x = x_good[i_c]
      y = y_good[i_c]
      
      #print np.min(np.sqrt((x-xlist)**2+(y-ylist)**2))
      if np.min(np.sqrt((x-xlist)**2+(y-ylist)**2)) < 2.*fwhm_beam_pix:
        #print 'prob x,y'
        result = False
        #print 'too close'
        #continue
      #if np.isnan(resid[x,y]):
        ## source is in a bad place on the map
        ##print 'prob nan'
        #result = False
        #continue
      if np.isnan(np.sum(resid[x-fwhm_beam_pix:x+fwhm_beam_pix,y-fwhm_beam_pix:y+fwhm_beam_pix])):
        # source is in a bad place on the map - somewhere a nan in a box fwhm_beam_pix  X fwhm_beam_pix 
        #print 'prob nan'
        result = False
        #print 'near nan'
        #continue
    #print  ' * ind %f mus' %((time.time()-T1)*1e6)
          
    # remove nearby pixels from good list for source coords
    #good_ind = np.where(np.sqrt((x-x_good)**2+(y-y_good)**2) > 2.*fwhm_beam_pix )
    #good_ind = np.where(np.abs(x-x_good) > 2.*fwhm_beam_pix )
    #x_good = x_good[good_ind]
    #y_good = y_good[good_ind]
    #good_ind = np.where(np.abs(y-y_good) > 2.*fwhm_beam_pix )
    #x_good = x_good[good_ind]
    #y_good = y_good[good_ind]
      
      
      # draw flux
    #print  ' * ind %f mus' %((time.time()-T1)*1e6)
    t1 = time.time()
    #flux = random_flux(Fmin=Fmin, Fmax=Fmax)
    #print 'get flux'
    flux = random_flux_powerlaw(alpha=alpha-1, Fmin=Fmin, Fmax=Fmax)
      
    pixcrd = np.array([[y,x,0,0]], np.float_)

    world = mapwcs.wcs_pix2sky(pixcrd, 0)
    
    t2 = time.time()
    #print  ' * %f mus' %((time.time()-T1)*1e6)
          
      
    #print  ' * %f mus' %((time.time()-T1)*1e6)
    if point:
        t1 = time.time()
        #print  ' ps %i of %i' %(i,Nsrc)
        sx = sig_beam_maj
        sy = sig_beam_min
        ba = 0.
        #pa = np.random.uniform(0.,180.)
        #pa = 90
        pa = beam_pa
        peak = flux #/(np.pi*sx*sy)   # Jy/beam in pix
        
        #sx = np.sqrt(2.)*np.sqrt(sig_beam**2.)
        #sy = np.sqrt(2.)*np.sqrt(sig_beam**2.)
        #pa = 0.
        #t1 = time.time()
        #model = add_gaussian_point(model,peak,x,y,sig_beam,gaus_point)
        #xmax,ymax =  data.shape
        x1 = np.floor(x) - Psize
        x2 = np.ceil(x) + Psize
        y1 = np.floor(y) - Psize
        y2 = np.ceil(y) + Psize
        
        try:
            model[x1:x2,y1:y2] += peak*gaus_point
            #t2 = time.time()
            #print 'p', t2-t1
            pnt[i] = 1
            t2 = time.time()
            #print  ' pd %i of %i %f ms' %(i,Nsrc,(t2-t1)*1e3)
        except:
            print x1,x2,y1,y2
            
    else:
        t1 = time.time()
        #print  ' es %i of %i' %(i,Nsrc)
        # more sources slightly resolved out to 2x beam
        sx = np.random.triangular(sig_beam,sig_beam, 4.*sig_beam)
        ba = np.random.triangular(0.4,1.,1.)  #rounder sources more common and no extreme
        sy = sx*ba
        if sy < sig_beam_min:
          sy = sig_beam_min
      
        pa = np.random.uniform(0.,180.)
        #pa = beam_bpa
        peak = flux /(sx*sy/(sig_beam_maj*sig_beam_min))   # Jy/beam in pix
        # convolve axes with beam
        
        #sx = np.sqrt(2*sx**2.) #+ sig_beam**2.)
        #sy = np.sqrt(2*sy**2.) #+ sig_beam**2.)
        
        #t1 = time.time()
        #model = add_gaussian(model,peak,x,y,sx,sy,pa)
        #add_gaussian(data,Speak,x0,y0,dx,dy,bpa):
        #print flux,x0,y0,dx,dy,bpa
        #Speak = flux/(2.*np.pi*dx*dy)
        #result = True
        #ymin = 0
        C = 6
        if sx/sy < 0.2:
            C = 8
        #xmax,ymax =  data.shape
        size = np.ceil(max(C*sx, C*sy))
        x1 = np.floor(x) - size
        x2 = np.ceil(x) + size
        y1 = np.floor(y) - size
        y2 = np.ceil(y) + size
        #print x1,x2
        #print y1,y2
        #for xi in range( x1, x2 ):
            #for yi in range( y1, y2 ):
            ##print xi,yi
            #data[xi,yi] += Speak*twodgaussian(xi,yi,1,x0,y0,dx,dy,bpa)
            #if resid[xi,y1] == 0.:
                #result = True
                #return data, result
            
        xg = np.arange(0, 2*size, 1, float)
        yg = xg[:,np.newaxis]
        rota = np.pi/180. * (pa-90.)
        x0 = size
        y0 = size
        x0p = x0 * np.cos(rota) - y0 * np.sin(rota)
        y0p = x0 * np.sin(rota) + y0 * np.cos(rota)
                    
        #print xg.shape, yg.shape, rota
        xp = xg * np.cos(rota) - yg * np.sin(rota)
        yp = xg * np.sin(rota) + yg * np.cos(rota)
        #yp = yp[:,np.newaxis]
        g = np.exp( -( (((x0p-xp)**2.)/(2.*sx**2.)) + (((y0p-yp)**2.)/(2.*sy**2.)) ) )
        
        
        #
        try:
            model[x1:x2,y1:y2] += peak*g
        except:
            print "ERROR!!"
            print size, x1,x2,y1,y2, model.shape, model[x1:x2,y1:y2].shape, g.shape, Speak
        #t2 = time.time()
        #print 'e', t2-t1
        n_ext2 += 1
        pnt[i] = 0
        t2 = time.time()
        #print  ' ed %i of %i %f ms' %(i,Nsrc,(t2-t1)*1e3)
    #print  ' * %f mus' %((time.time()-T1)*1e6)
    
    T2 = time.time()
    #print 'T = %f ms' %((T2-T1)*1e3)
        
    xlist[i] = x
    ylist[i] = y
    Speak[i] = peak
    Flux[i] = flux #/(np.pi*fwhm_beam**2. ))
    #bmaj[i] = sx*sig_beam*ps  #in arcsec
    #bmin[i] = sy*sig_beam*ps  #in arcsec
    bmaj[i] = sx*fwhm_f*ps  #in arcsec
    bmin[i] = sy*fwhm_f*ps  #in arcsec
    bpa[i] = pa
    ra[i] = world[0][0] #-ps/3600.
    dec[i] = world[0][1] #-ps/3600.
    #print x,y,ra[i], dec[i], world
  # end generate sources #
    
  print '%i extended sources (%.3f per cent) ' %(n_ext2, 100.*n_ext2/Nsrc)
 
  # add noise to map
  print 'adding noise'
  modeln = model+ resid
  
  # write to fits
  modelfits = model_map_name
  if os.path.isfile( modelfits ):
    os.system(  'rm '+modelfits )
  #inhead['NAXIS1'] = int(nx)
  #inhead['NAXIS2'] = int(ny)
  #modeln = np.array([[modeln]])
  #hdu = pf.PrimaryHDU(modeln, inhead)
  #hdu.writeto(modelfits)
  #tt=t[:,:,np.newaxis, np.newaxis]
  pf.writeto(modelfits, modeln[np.newaxis,np.newaxis,:,:], sourcehead)
  
  # write input catalogue
  catfile = file(model_cat+'.tempascii','w')
  catfile.write('# Source_id ra dec Total_flux Peak_flux Bmaj Bmin Bpa\n')
  #print len(Speak)
  for i in range(len(Speak)):
    catfile.write('%i %.6f %.6f %f %f %f %f %f\n' %(i,ra[i],dec[i],Flux[i],Speak[i],bmaj[i],bmin[i],bpa[i]))
  catfile.close()
  
  # Region file format: DS9 version 4.0
  catfile = file(model_cat+'.reg','w')
  head = '''global color=yellow font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source
fk5
'''
  catfile.write(head)
  for i in range(len(Speak)):
    if pnt[i] == 1:
      catfile.write('ellipse(%f,%f,%f",%f",%f)# text={%i} color={yellow}\n' %(ra[i],dec[i],bmaj[i],bmin[i],bpa[i],i))
    else:
      catfile.write('ellipse(%f,%f,%f",%f",%f)# text={%i} color={blue}\n' %(ra[i],dec[i],bmaj[i],bmin[i],bpa[i]-90.,i))
  catfile.close()

  
  # save as fits table
  os.system("stilts tcopy in=%s out=%s ifmt=ascii ofmt=fits" %(model_cat+'.tempascii',model_cat))
  #os.system("rm -rf %s" %(model_cat+'.tempascii'))
  
  return Speak, Flux
    

def get_flux_slope(fluxes,rms):
    
    n,bins = np.histogram(np.log10(fluxes),bins=25,range=(np.log10(5*rms), max(np.log10(fluxes))))
    bin_c = (bins[1:]+bins[0:-1])/2
    bin_c = bin_c[n!=0]
    n = n[n!=0]
    dbin = np.average(bins[1:]-bins[0:-1])
    nperbin = n/dbin
    nperbin = np.log10(nperbin)
    p = np.polyfit(bin_c, nperbin, 1.)
    x = np.linspace(bin_c.min(), bin_c.max(), 100)
    y = np.polyval(p,x)
    #plot(bin_c, nperbin, 'b')
    #plot(x,y,'r')
    
    slope = p[0]
    
    return slope  

def get_areas(catfile,cat_areafile, rmsfits = 'RMS.FITS'):
  print 'calculating source detection areas'
  rmsdat = pf.getdata(rmsfits)
  rmsdat = np.ma.masked_where(np.isnan(rmsdat), rmsdat).compressed()
  rmshead = pf.getheader(rmsfits)
  ps = rmshead.get('CDELT2')
  cat = read_match_cat(catfile)
  cat_areas = file(cat_areafile,'w')
  cat_areas.write('Source_id,AreaP,AreaS,AreaPo,AreaSo\n')
  names = cat[0]
  flux = []
  sources = []
  pi = names.index('peak_flux')
  pi2 = names.index('peak_flux', names.index('peak_flux')+1)
  si = names.index('source_id')
  rmsmin = rmsdat.min()
  for i, entry in enumerate(cat[1:]):
    source = entry[si]
    fpi = entry[pi]
    fpo = entry[pi2]
    #if fpi < 5.*rmsmin:
      #pixels = 0.
      #area = 0.
    #else:
    if 1:
      rmslim = fpi
      pixels = np.sum(rmsdat < rmslim)
      area = pixels*ps*ps
    if fpo < 0:
      pixelso = -1
      areao = -1
    #elif fpo < 5.*rmsmin:
      #pixelso = 0.
      #areao = 0.
    else:
      rmslimo = fpo
      pixelso = np.sum(rmsdat < rmslimo)
      areao = pixelso*ps*ps
    cat_areas.write('%f,%f,%f,%f,%f\n' %(source,pixels,area,pixelso,areao))
  cat_areas.close()
  print 'done'
  return
  
def clip_mean_sig(x,w,sigclip=3.,N=10.):
    av = np.average(x,weights=w)
    #s = np.dot( w, (x-av)**2./w.sum() )
    s = np.std(x)
    for i in np.arange(N):
      mask = abs(x-av) > sigclip*s
      x = np.ma.masked_where(mask,x).compressed()
      w = np.ma.masked_where(mask,w).compressed()
      av = np.average(x,weights=w)
      #s = np.dot( w, (x-av)**2./w.sum() )
      s = np.std(x)
    return av,s
def clip_mean_sig_nw(x,sigclip=3.,N=10.):
    av = np.average(x)
    #s = np.dot( w, (x-av)**2./w.sum() )
    s = np.std(x)
    for i in np.arange(N):
      mask = abs(x-av) > sigclip*s
      x = np.ma.masked_where(mask,x).compressed()
      av = np.average(x)
      #s = np.dot( w, (x-av)**2./w.sum() )
      s = np.std(x)
    return av,s
    


def file_exists( name ):
  return os.path.isfile(name)
  
def find_sources(fitsdir, fitsname, catname, clobber=False, thresh_pix=5.0, thresh_isl=3.0, rms_box=(80,40), output_files=False):  
  # run PyBDSM
  if (file_exists( catname ) and clobber) or (not file_exists( catname )):
    pwd = os.getcwd()
    os.chdir(fitsdir)
    bdsm_img = bdsm.process_image( fitsname, adaptive_rms_box=False, rms_map = True, rms_box = rms_box, advanced_opts = True, group_by_isl = True,  thresh_pix = thresh_pix, thresh_isl = thresh_isl ) #, rms_map = True, rms_box = (100,25))
    
    bdsm_img.write_catalog( format = 'fits', outfile = catname.replace('.fits','.gaul.fits'), clobber=True, catalog_type = 'gaul', incl_empty=False )
    bdsm_img.write_catalog( format = 'ds9', outfile = catname.replace('.fits','.gaul.ds9'), clobber=True, catalog_type = 'gaul', incl_empty=False )
    bdsm_img.write_catalog( format = 'fits', outfile = catname, clobber=True, catalog_type = 'srl', incl_empty=False )
    bdsm_img.write_catalog( format = 'ds9', outfile = catname.replace('.fits','.ds9'), clobber=True, catalog_type = 'srl', incl_empty=False )
    f = file(catname.replace('.fits','.ds9'),'r')
    f2 = file(catname.replace('.fits','.reg'),'w')
    for line in f:
        f2.write(line.replace(fitsname.replace('.fits',''),''))
    f.close()
    f2.close()
    
    if output_files:
        rmsfits = fitsname.replace('.FITS','.RMS.FITS')
        rmsfits = fitsname.replace('.fits','.rms.fits')
        bdsm_img.export_image(outfile = rmsfits, img_type='rms', clobber=True)
        
        resfits = fitsname.replace('.FITS','.RES.FITS')
        resfits = fitsname.replace('.fits','.res.fits')
        bdsm_img.export_image(outfile = resfits, img_type='gaus_resid', clobber=True)
    
    
    os.chdir( pwd )
  return
    
def find_sources_cmd(fitsdir,findcmd):  
    
  pwd = os.getcwd()
  os.chdir(fitsdir)
  
  print pwd
  print fitsdir
  
  print findcmd
  os.system(findcmd)
  
  os.chdir( pwd )
  
  ## run PyBDSM
  #if (file_exists( catname ) and clobber) or (not file_exists( catname )):
    #pwd = os.getcwd()
    #os.chdir(fitsdir)
    #bdsm_img = bdsm.process_image( fitsname, adaptive_rms_box=False, rms_map = True, rms_box = rms_box, advanced_opts = True, group_by_isl = True,  thresh_pix = thresh_pix, thresh_isl = thresh_isl ) #, rms_map = True, rms_box = (100,25))
    
    #bdsm_img.write_catalog( format = 'fits', outfile = catname.replace('.fits','.gaul.fits'), clobber=True, catalog_type = 'gaul', incl_empty=False )
    #bdsm_img.write_catalog( format = 'ds9', outfile = catname.replace('.fits','.gaul.ds9'), clobber=True, catalog_type = 'gaul', incl_empty=False )
    #bdsm_img.write_catalog( format = 'fits', outfile = catname, clobber=True, catalog_type = 'srl', incl_empty=False )
    #bdsm_img.write_catalog( format = 'ds9', outfile = catname.replace('.fits','.ds9'), clobber=True, catalog_type = 'srl', incl_empty=False )
    #f = file(catname.replace('.fits','.ds9'),'r')
    #f2 = file(catname.replace('.fits','.reg'),'w')
    #for line in f:
        #f2.write(line.replace(fitsname.replace('.fits',''),''))
    #f.close()
    #f2.close()
    
    #if output_files:
        #rmsfits = fitsname.replace('.FITS','.RMS.FITS')
        #rmsfits = fitsname.replace('.fits','.rms.fits')
        #bdsm_img.export_image(outfile = rmsfits, img_type='rms', clobber=True)
        
        #resfits = fitsname.replace('.FITS','.RES.FITS')
        #resfits = fitsname.replace('.fits','.res.fits')
        #bdsm_img.export_image(outfile = resfits, img_type='gaus_resid', clobber=True)
    
    
    #os.chdir( pwd )
  return
    

def mask_bad_resid(newresidmap, residmap, rmsmap, thresh=0.025):
    print 'masking residual map'
    inhead = pf.getheader(residmap)
    resid = pf.getdata(residmap)[0,0,:,:]
    resid2 = resid.transpose()
    rms = pf.getdata(rmsmap)[0,0,:,:]
    #mask = (rms > thresh)
    mask_val = np.where((rms > thresh))
    mask_x = mask_val[0]
    mask_y = mask_val[1]
    nx,ny = resid.shape
    mask_resid = resid
    for i in range(len(mask_x)):
        xi = mask_x[i]
        yi = mask_y[i]
        if not np.isnan(resid2[xi,yi]):
            mask_resid[xi,yi] = resid2[xi,yi]
        elif not np.isnan(resid2[xi,xi]):
            mask_resid[xi,yi] = resid2[xi,xi]
        elif not np.isnan(resid2[yi,yi]):
            mask_resid[xi,yi] = resid2[yi,yi]
    #for xi in range(nx):
        #for yi in range(ny):
            #if rms[xi,yi] > thresh:
                #mask_resid[xi,yi] = resid2[xi,yi]
    mask_val = np.where((resid > 3.*rms.max()))
    mask_x = mask_val[0]
    mask_y = mask_val[1]
    nx,ny = resid.shape
    mask_resid = resid
    for i in range(len(mask_x)):
        xi = mask_x[i]
        yi = mask_y[i]
        if not np.isnan(resid2[xi,yi]):
            mask_resid[xi,yi] = resid2[xi,yi]
        elif not np.isnan(resid2[xi,xi]):
            mask_resid[xi,yi] = resid2[xi,xi]
        elif not np.isnan(resid2[yi,yi]):
            mask_resid[xi,yi] = resid2[yi,yi]
        
    hdu = pf.PrimaryHDU(mask_resid, inhead)
    hdu.writeto(newresidmap)
    print 'done'
    return


def snr_env2(p, x):
    return p[1]/x**p[2] + p[0]
def snr_env2_res(p, x, y):
    return y - snr_env2(p,x)

def flux_ratio_env_res (p,x,y):
  res = y - 1. - (p[0]**2. + (p[1]/x)**2.)**0.5
  return res
def flux_ratio_env (p,x):
  f = 1. + (p[0]**2. + (p[1]/x)**2.)**0.5
  return f


def time_smearing(resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # output and input units are equal
    # Condition that delta_Theta * delta_T << 1.37E4 * Resolution
    # where delta_Theta is the angular radius of the image and delta_T is the time smoothing.

    # Given resolution and field of view (delta_theta) this returns a condition on the time smearing at the correlator
    delta_T = 1.37E4*resolution/delta_Theta

    print 'Time averaging should be less than %s'%delta_T


    # Time average smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0

    print 'At radius %s a source will only have %s percent of its flux if data smoothed in time to %s'%(delta_Theta,Reduction,delta_T)
    
    return delta_T

def time_smearing2(delta_T,resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # output and input units are equal
    # Condition that delta_Theta * delta_T << 1.37E4 * Resolution
    # where delta_Theta is the angular radius of the image and delta_T is the time smoothing.
    # Time average smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B

    # Given resolution and field of view (delta_theta) this returns a condition on the time smearing at the correlator
    Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0

    return Reduction
    
def bandwidth_smearing(freq,resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal 
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Bandwidth smearing amplitude loss in http://adsabs.harvard.edu/full/1989ASPC....6..247B
    
    # Given resolution, freq and offset this gives the condition for the delta_freq

    delta_freq = freq*resolution/delta_Theta

    print 'Bandwidth averaging should be much less than %s'%delta_freq

    
    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    print 'At radius %s a source will only have %s percent of its flux if data smoothed in freq to %s'%(delta_Theta,Reduction,delta_freq)

    return delta_freq

def bandwidth_smearing2(delta_freq,freq,resolution,delta_Theta):
    #http://www.cv.nrao.edu/course/astr534/Interferometers1.html
    # Output and input units are equal 
    # Condition that delta_Theta * delta_freq << Freq * resolution
    # where delta_Theta is the offset from the pointing centre and delta_freq is the bandwidth smoothing.

    # Given resolution, freq and offset this gives the condition for the delta_freq

    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    return Reduction



def get_peak_total_env_frac(OutCatList, SOURCElist, SNRrange, title='', savename='', ylim=[0.51,25.], logy=True):
    import astropy.coordinates as ac
    
    outpeak_all =  np.array([])
    outflux_all = np.array([])
    rms_all=np.array([])
    loss_g_all=np.array([])
    ra_all=np.array([])
    dec_all=np.array([])
    
    f2, ax2 = paper_single_ax()
    for OutCat in OutCatList:
        
        cat_out = load_fits(OutCat)
        print OutCat, ": ", len(cat_out), " sources"
        
        cat_out = cat_out[cat_out.Total_flux > 0]
        print OutCat, ": ", len(cat_out), " sources"
        
        outpeak = cat_out.Peak_flux.copy()
        outflux = cat_out.Total_flux.copy()
        rms = cat_out.Isl_rms.copy()
              

        ra_all = np.hstack((ra_all, cat_out['RA']))
        dec_all = np.hstack((dec_all, cat_out['DEC']))
        
        outpeak_all = np.hstack((outpeak_all, outpeak))
        outflux_all = np.hstack((outflux_all, outflux))
        rms_all = np.hstack((rms_all, rms))
        
        #plot(cat_out['RA'], cat_out['DEC'], 'ko',alpha=0.1)
        
        loss2 = 0.833
        loss2 = 1.
        RA0 = 218.
        DEC0 = 34.5
        rad_g = angular_separation(cat_out['RA'], cat_out['DEC'],RA0,DEC0)
        bloss_g = bandwidth_smearing2(97.5e3,150.e6,7.,rad_g*3600.)
        tloss_g = time_smearing2(8, 7.,rad_g*3600)
        loss_g = bloss_g*tloss_g*loss2
    
        loss_g_all = np.hstack((loss_g_all, loss_g))
    
    ## get rid of the artefacts common to all maps
    catc=ac.SkyCoord(ra_all, dec_all, unit='deg')
    t,ang,d = ac.match_coordinates_sky(catc,catc,nthneighbor=2)
    goodind = ang.arcsec > 0.001
    outpeak_all = outpeak_all[goodind]
    outflux_all = outflux_all[goodind]
    rms_all = rms_all[goodind]
    
    
    SNR = outpeak_all/rms_all
    Srat = outflux_all/(outpeak_all)
    
    print len(Srat)," total sources"

    
    import matplotlib.pyplot as plt
    
    print '## peak/total ratio ##'
    #f2 = plt.figure(figsize = (7.0,5.0))
    #ax2 = plt.subplot(111)
    f2, ax2 = paper_single_ax()
    ax2.semilogx(SNR, Srat, '.', color='gray',alpha=0.1)

    env = np.ones(len(SNRrange))*1.02
    #nsrci = np.sum(detect)
    for k in range(len(SNRrange)-1):
        SiSp = np.ma.masked_where(1-((SNR>SNRrange[k])*(SNR<SNRrange[k+1])),Srat).compressed()
        #print '#',SNRrange[k], SNRrange[k+1], len(SiSp)
        if len(SiSp) > 1:
            SiSpi = 1.0
            #print i,SiSpi, np.sum(SiSp<SiSpi), np.sum(SiSp<SiSpi)/len(SiSp)
            while float(np.sum(SiSp<SiSpi))/len(SiSp) < 0.95:
                #print i,SiSpi, np.sum(SiSp<SiSpi), np.sum(SiSp<SiSpi)/len(SiSp)
                #print k, SNRrange[k], len(SiSp),np.sum(SiSp<SiSpi)
                SiSpi+= 0.01
            env[k] = SiSpi
    env[ k+1 ] = env[ k ]
    #rms_r[k] = rms

    #snr_S.append(SNRrange)
    #fraction_S.append(env)
    #ax2.semilogx(SNRrange[:-4], env[:-4], 'ko')
    ax2.semilogx(SNRrange, env, 'ko', label=r'95\% envelope')
    #ax2.vlines(SNRrange,0,6,'k')
    
    #p0 = [0.04,2.5]
    p0 = [1., 150., 0.5]
    print 'fitting 95 percent envelope'
    from scipy.optimize import leastsq
    #plsq,cov,info,mesg,success = leastsq(flux_ratio_env_res, p0, args=(SNRrange, env), full_output=True)
    plsq,cov,info,mesg,success = leastsq(snr_env2_res, p0, args=(SNRrange, env), full_output=True)
    if success==1:
        # calculate final chi square
        chisq=sum(info["fvec"]*info["fvec"])
        dof=len(SNRrange)-len(p0)
        # chisq, sqrt(chisq/dof) agrees with gnuplot
        print "Converged with chi squared ",chisq
        print "degrees of freedom, dof ", dof
        print "RMS of residuals (i.e. sqrt(chisq/dof)) ", np.sqrt(chisq/dof)
        print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
        # uncertainties are calculated as per gnuplot, "fixing" the result
        # for non unit values of the reduced chisq.
        # values at min match gnuplot
        print "Fitted parameters at minimum, with 68% C.I.:"
        for i,pmin in enumerate(plsq):
            try:
                print "%2i %12f +/- %10f"%(i,pmin,np.sqrt(cov[i,i])*np.sqrt(chisq/dof))
            except:
                print i,pmin,np.sqrt(cov)*np.sqrt(chisq/dof)
        print "Correlation matrix"
        # correlation matrix close to gnuplot
        #for i in range(len(plsq)):
            #print "%10i"%i,
            #for j in range(i+1):
                #try:
                    #print "%10f"%(cov[i,j]/np.sqrt(cov[i,i]*cov[j,j]),),
                #except:
                    #print 
            #print
    else: 
        print "Not converged", mesg
    #-----------------------------------------------


    print '%f, %f' %(plsq[0], plsq[1])
    plt.xlim(2e0,1.2e3)
    SNRrange_env = np.logspace(0.,3.,100)
    #envelope = flux_ratio_env(plsq[0], SNRrange_env )
    #envelope = 1. + (plsq[0]**2. + (plsq[1]/SNRrange_env)**2.)**0.5
    envelope = snr_env2(plsq, SNRrange_env)
    #envelope2 = 1. - (plsq[0]**2. + (plsq[1]/SNRrange_env)**2.)**0.5
    plt.semilogx(SNRrange_env,envelope,'r', label='Fit')
    #plt.semilogx(SNRrange_env,envelope2,'r:')

    
    ax2.set_ylabel('$S_{\mathrm{int}}/ S_{\mathrm{peak}}$')
    ax2.set_xlabel('$S_{\mathrm{peak}} / \sigma_L $ ')
    ax2.legend()
    ax2.semilogy()
    ax2.set_ylim(ylim[0], ylim[1])
    x1,x2 = ax2.get_xlim()
    ax2.set_xlim(4.5, x2)
    ax2.minorticks_on()
    plt.subplots_adjust(bottom=0.15,right=0.95)
    #f2.savefig('%s_peak_total_pointsim.png' %(name))
    ax2.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
    ax2.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    if logy:
        ax2.yaxis.set_major_locator(LogLocator(subs=[1,2,5]))
        ax2.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    else:
        ax2.yaxis.set_major_locator(MultipleLocator(1.0))
    fig_save_many(f2, '%s_peak_total_pointsim' %(savename),dpi=600)
    
    point_frac = plot_peak_total( savename, SOURCElist, plsq )
    
            
    return point_frac,  plsq



def get_size_dist(OutCatList, SOURCElist, SNRrange, title='', savename='', ylim=[0.51,25.], logy=True):
    import astropy.coordinates as ac
    
    plot_major( savename, SOURCElist )
    
    outpeak_all =  np.array([])
    outflux_all = np.array([])
    rms_all=np.array([])
    loss_g_all=np.array([])
    ra_all=np.array([])
    dec_all=np.array([])
    size_all=np.array([])
    
    f2, ax2 = paper_single_ax()
    for OutCat in OutCatList:
        
        cat_out = load_fits(OutCat)
        print OutCat, ": ", len(cat_out), " sources"
        
        cat_out = cat_out[cat_out.Total_flux > 0]
        print OutCat, ": ", len(cat_out), " sources"
        
        outpeak = cat_out.Peak_flux.copy()
        outflux = cat_out.Total_flux.copy()
        rms = cat_out.Isl_rms.copy()
        size = cat_out.Maj.copy()
        size = size *3600.  # in arcsec
              

        ra_all = np.hstack((ra_all, cat_out['RA']))
        dec_all = np.hstack((dec_all, cat_out['DEC']))
        
        outpeak_all = np.hstack((outpeak_all, outpeak))
        outflux_all = np.hstack((outflux_all, outflux))
        rms_all = np.hstack((rms_all, rms))
        size_all = np.hstack((size_all, size))
        
        #plot(cat_out['RA'], cat_out['DEC'], 'ko',alpha=0.1)
        
        loss2 = 0.833
        loss2 = 1.
        RA0 = 218.
        DEC0 = 34.5
        rad_g = angular_separation(cat_out['RA'], cat_out['DEC'],RA0,DEC0)
        bloss_g = bandwidth_smearing2(97.5e3,150.e6,7.,rad_g*3600.)
        tloss_g = time_smearing2(8, 7.,rad_g*3600)
        loss_g = bloss_g*tloss_g*loss2
    
        loss_g_all = np.hstack((loss_g_all, loss_g))
    
    ## get rid of the artefacts common to all maps
    catc=ac.SkyCoord(ra_all, dec_all, unit='deg')
    t,ang,d = ac.match_coordinates_sky(catc,catc,nthneighbor=2)
    goodind = ang.arcsec > 0.001
    outpeak_all = outpeak_all[goodind]
    outflux_all = outflux_all[goodind]
    rms_all = rms_all[goodind]
    size_all = size_all[goodind]
    
    
    SNR = outpeak_all/rms_all
    Srat = size_all
    
    print len(Srat)," total sources"

    
    import matplotlib.pyplot as plt
    
    print '## peak/total ratio ##'
    #f2 = plt.figure(figsize = (7.0,5.0))
    #ax2 = plt.subplot(111)
    f2, ax2 = paper_single_ax()
    ax2.semilogx(SNR, Srat, '.', color='gray',alpha=0.1)

    plt.xlim(2e0,1.2e3)

    
    ax2.set_ylabel('$a$ [arcsec]')
    ax2.set_xlabel('$S_{\mathrm{peak}} / \sigma_L $ ')
    ax2.legend()
    #ax2.semilogy()
    ax2.set_ylim(ylim[0], ylim[1])
    x1,x2 = ax2.get_xlim()
    ax2.set_xlim(4.5, x2)
    ax2.minorticks_on()
    plt.subplots_adjust(bottom=0.15,right=0.95)
    #f2.savefig('%s_peak_total_pointsim.png' %(name))
    ax2.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
    ax2.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    if logy:
        ax2.yaxis.set_major_locator(LogLocator(subs=[1,2,5]))
        ax2.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    else:
        ax2.yaxis.set_major_locator(MultipleLocator(1.0))
        
    catname = SOURCElist
    cat = load_fits(catname)
    
    flux = cat.Total_flux
    eflux = cat.E_Total_flux
    peak = cat.Peak_flux
    epeak = cat.E_Peak_flux
    rms = cat.Isl_rms
    major = cat.Maj*3600 # in arcsec
    emajor = cat.E_Maj*3600  #in arcsec
    
    
    #f2, ax2 = paper_single_ax()
    x=peak/rms
    y=major
    plt.semilogx(x,y,'k.',alpha=0.2, mec='none')
    sy = emajor
    plt.xlim(2e0,1.2e3)
    SNRrange = np.logspace(0.,3.,100)
    
    envy=flux/peak
    
    plt.ylabel(r'$a$ [arcsec]')
    plt.xlabel(r'$S_{\mathrm{peak}}/\sigma_L$')
    
    plt.minorticks_on()
    plt.subplots_adjust(bottom=0.15,right=0.95)
    # force ylimits
    #plt.ylim(2.,10.)
    print '** %i sources  Si/Sp>6' %(sum(y>6))
    #f2.savefig('%s_peak_total.png' %(name))
    #f2.savefig('%s_peak_total.eps' %(name),dpi=600)
    ax2.yaxis.set_major_locator(MultipleLocator(10.0))
    ax2.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
    ax2.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    #fig_save_many(f2, '%s_major' %(name),dpi=600)
        
    fig_save_many(f2, '%s_size_extsim' %(savename),dpi=600)
    
            
    return 

def plot_peak_total( name, catname , p , logy=True):
  import matplotlib.pyplot as plt
  
  cat = load_fits(catname)
  
  flux = cat.Total_flux
  eflux = cat.E_Total_flux
  peak = cat.Peak_flux
  epeak = cat.E_Peak_flux
  rms = cat.Isl_rms
  
  
        
  loss2 = 0.833
  loss2 = 1
  RA0 = 218.
  DEC0 = 34.5
  rad_g = angular_separation(cat['RA'], cat['DEC'],RA0,DEC0)
  bloss_g = bandwidth_smearing2(97.5e3,150.e6,7.,rad_g*3600.)
  tloss_g = time_smearing2(8, 7.,rad_g*3600)
  loss_g = bloss_g*tloss_g*loss2

  
  f2, ax2 = paper_single_ax()
  x=peak/rms
  y=flux/(peak/loss_g)
  plt.semilogx(x,y,'k.', mec='none', alpha=0.2)
  sy = np.sqrt( (epeak*y/peak)**2. + (eflux*y/flux)**2. )
  plt.xlim(2e0,1.2e3)
  SNRrange = np.logspace(0.,3.,100)
  envelope = snr_env2(p, SNRrange)
  plt.semilogx(SNRrange,envelope,'r', ls='dashed', label='Fitted envelope')
  
  count_resolved = 0
  count_unresolved = 0
  for i in range(len(peak)):
    if y[i] > snr_env2(p, peak[i]/rms[i]):
    #if (y[i]-abs(sy[i])) > (1+ ((p[0])**2.+(p[1]/(peak[i]/rms[i]))**2.)**0.5):
      count_resolved += 1
    else:
      count_unresolved += 1
  print '# resolved = %i' %(count_resolved)
  print '# unresolved = %i' %(count_unresolved)
  
  point_frac = 1.-1.*count_resolved/(count_resolved+count_unresolved)
  
  print '# point fraction = %.3f' %(point_frac)
  
  #plt.title(' Source Counts ')
  plt.ylabel(r'$S_{\mathrm{int}} /S_{\mathrm{peak}}$')
  plt.xlabel(r'$S_{\mathrm{peak}}/\sigma_L$')
  
  plt.minorticks_on()
  plt.subplots_adjust(bottom=0.15,right=0.95)
  # force ylimits
  #plt.ylim(0.,6.)
  plt.semilogy()
  plt.ylim(0.51,25.)
  print '** %i sources  Si/Sp>6' %(sum(y>6))
  ax2.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
  if logy:
      ax2.yaxis.set_major_locator(LogLocator(subs=[1,2,5]))
      ax2.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
  else:
      ax2.yaxis.set_major_locator(MultipleLocator(1.0))
  ax2.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
  ax2.legend()
  x1,x2 = ax2.get_xlim()
  ax2.set_xlim(4.5, x2)
  fig_save_many(f2, '%s_peak_total' %(name),dpi=600)
  return point_frac

def plot_major( name, catname ):
  import matplotlib.pyplot as plt
  
  cat = load_fits(catname)
  
  flux = cat.Total_flux
  eflux = cat.E_Total_flux
  peak = cat.Peak_flux
  epeak = cat.E_Peak_flux
  rms = cat.Isl_rms
  major = cat.Maj*3600 # in arcsec
  emajor = cat.E_Maj*3600  #in arcsec
  
  
  f2, ax2 = paper_single_ax()
  x=peak/rms
  y=major
  plt.semilogx(x,y,'k.',alpha=0.2, mec='none')
  sy = emajor
  plt.xlim(2e0,1.2e3)
  SNRrange = np.logspace(0.,3.,100)
  
  envy=flux/peak
  
  plt.ylabel(r'$a$ [arcsec]')
  plt.xlabel(r'$S_{\mathrm{peak}}/\sigma_L$')
  
  plt.minorticks_on()
  plt.subplots_adjust(bottom=0.15,right=0.95)
  # force ylimits
  #plt.ylim(2.,10.)
  print '** %i sources  Si/Sp>6' %(sum(y>6))
  #f2.savefig('%s_peak_total.png' %(name))
  #f2.savefig('%s_peak_total.eps' %(name),dpi=600)
  ax2.yaxis.set_major_locator(MultipleLocator(10.0))
  ax2.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
  ax2.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
  fig_save_many(f2, '%s_major' %(name),dpi=600)
  return 

def plot_bright_size_distrib( name, catname , snrlim=20 ):
  import matplotlib.pyplot as plt
 
  cat = load_fits(catname)
  
  flux = cat.Total_flux
  eflux = cat.E_Total_flux
  peak = cat.Peak_flux
  epeak = cat.E_Peak_flux
  rms = cat.Isl_rms
  maj = cat.Maj*3600  # in arcsec
  
  
  f2, ax2 = paper_single_ax()
  #plt.figure(figsize=(7.0,5.0))
  x=peak/rms
  y=flux/peak

  ax2.hist(maj, bins=20, label='all', histtype='step')
  ax2.hist(maj[x>snrlim], label='selection SNR<'+r'{lim}'.format(lim=snrlim), bins=20, histtype='step')
  ax2.legend()

  plt.title(' Major axis for sources')
  plt.ylabel(r'count')
  plt.xlabel(r'Maj [arcsec]')
  
  plt.minorticks_on()
  plt.subplots_adjust(bottom=0.15,right=0.95)

  fig_save_many(f2, '%s_size distrib' %(name),dpi=600)
  return 

def plot_bright_size_distrib( name, catname , snrlim=20 ):
  import matplotlib.pyplot as plt
 
  cat = load_fits(catname)
  
  flux = cat.Total_flux
  eflux = cat.E_Total_flux
  peak = cat.Peak_flux
  epeak = cat.E_Peak_flux
  rms = cat.Isl_rms
  maj = cat.Maj*3600  # in arcsec
  
  
  f2, ax2 = paper_single_ax()
  #plt.figure(figsize=(7.0,5.0))
  x=peak/rms
  y=flux/peak

  ax2.hist(maj, bins=20, label='all', histtype='step')
  ax2.hist(maj[x>snrlim], label='selection SNR<'+r'{lim}'.format(lim=snrlim), bins=20, histtype='step')
  ax2.legend()

  plt.title(' Major axis for sources')
  plt.ylabel(r'count')
  plt.xlabel(r'Maj [arcsec]')
  
  plt.minorticks_on()
  plt.subplots_adjust(bottom=0.15,right=0.95)

  fig_save_many(f2, '%s_size distrib' %(name),dpi=600)
  return 

# make fake images

#def make_model_maps(N, modname, Fmin, Fmax, residmap, Nsrc,  verbose=False):
    #if isinstance(N, int):
        #iterate = range(N)
    #elif isinstance(N, list):
        #if len(N) == 2:
            #iterate = range(N[0], N[1])
    #else:
        #raise ValueError 
    
    #maplist = []
    #catlist = []
    #for i in iterate:
        #if verbose:  print 'RUN %i' %(i)
        #model_map_name = '{modname:s}_{i1:i}.fits'.format(modname,i1=i)
        #model_map_cat_name = '{modname:s}_{i1:i}.modelcat.in.fits'.format(modname,i1=i)

        #exists = os.path.isfile(MODEL_DIR+model_map_name)
        #if (exists and clobber_model) or (not exists):
            
            ##0.1*FLUXMIN*1000., 0.75*FLUXMAX*1000.
            #speak,flux = make_model_map(MODEL_DIR+model_map_name, MODEL_DIR+model_map_cat_name, Fmin=Fmin, Fmax=Fmax, residmap=residmap, Nsrc=Nsrc)
            
        #maplist.append(model_map_name)
        #catlist.append(model_map_cat_name)
    #return maplist, catlist
        

#modelmaps, modelincats = make_model_maps([Nsim0,NsimN], modname=FIELDNAME, Fmin=0.1*FLUXMIN, Fmax=0.75*FLUXMAX, residmap = RESIDUALmap, Nsrc=NSRC,  verbose=False)


def fix_beam_in_header(infitsname):
    infits = pf.open(infitsname, mode='update')
    inhead = infits[0].header
    beam_maj = inhead.get('BMAJ')
    beam_min = inhead.get('BMIN')
    # or it is not there, check the history
    hist = inhead.get_history()
    if beam_maj is None:
        for i in range(len(hist)):
            t = hist[i]
            if 'BMAJ=' in t:
                tt = t.split()
                beam_maj = float(tt[tt.index('BMAJ=')+1])
                break
    if beam_min is None:
        for i in range(len(hist)):
            t = hist[i]
            if 'BMIN=' in t:
                tt = t.split()
                beam_min = float(tt[tt.index('BMIN=')+1])
                break
    inhead['BMAJ'] = beam_maj
    inhead['BMIN'] = beam_min
    infits.flush()
    return


def plot_input_output_cat(inputcat, outputcat, savename, bins=20):
    import matplotlib.pyplot as plt
    
    incat = load_fits(inputcat)
    outcat = load_fits(outputcat)
    
    f, ax = paper_double_mult_ax(nrows=2,ncols=3,setticks=False,sharey=True)
    ax1 = ax[0][0]
    ax2 = ax[0][1]
    ax3 = ax[0][2]
    
    n,b,p = ax1.hist(np.log10(incat.Total_flux*1000), bins=bins,label='Input',histtype='step')
    ax1.hist(np.log10(outcat.Total_flux*1000), bins=b, label='Output',histtype='step')
    ax1.set_xlabel('log Total Flux [mJy]')
    ax1.set_ylabel('count')
    
    n,b,p = ax2.hist(np.log10(incat.Peak_flux*1000), bins=bins,label='Input',histtype='step')
    ax2.hist(np.log10(outcat.Peak_flux*1000), bins=b, label='Output',histtype='step')
    ax2.set_xlabel('log Peak Flux [mJy]')
    
    n,b,p = ax3.hist(incat.Peak_flux/incat.Total_flux, bins=bins,label='Input',histtype='step')
    ax3.hist(outcat.Peak_flux/outcat.Total_flux, bins=b, label='Output',histtype='step')
    ax3.set_xlabel('Peak/Total Flux')
    ax3.legend()
    
    #f.savefig(savename+'_sourcefluxes.png')
    
    plt.close(f)
    
    #f, ax = paper_double_mult_ax(nrows=2,ncols=3,setticks=False,sharey=True)
    #f, (ax1, ax2, ax3) = plt.subplots(1,2,sharey=True)
    ax1 = ax[1][0]
    ax2 = ax[1][1]
    ax3 = ax[1][2]
    
    n,b,p = ax1.hist(incat.Bmaj, bins=bins ,label='Input',histtype='step')
    ax1.hist(outcat.Maj*3600, bins=bins, label='Output',histtype='step')
    ax1.set_xlabel('Bmaj [arcsec]')
    ax1.set_ylabel('count')
    
    n,b,p = ax2.hist(incat.Bmin, bins=bins ,label='Input',histtype='step')
    ax2.hist(outcat.Min*3600, bins=bins, label='Output',histtype='step')
    ax2.set_xlabel('Bmin [arcsec]')
    #ax2.set_ylabel('count')
    
    n,b,p = ax3.hist(incat.Bpa, bins=bins ,label='Input',histtype='step')
    ax3.hist(outcat.PA, bins=bins, label='Output',histtype='step')
    ax3.set_xlabel('Bpa')
    #ax3.set_ylabel('count')
    ax3.legend()
    
    f.savefig(savename)
    
    plt.close(f)
    
    return


class SimMap:
    
    def __init__(self, inputmap, inputcat, residualmap, pbmap, N, Fmin, Fmax, dirname, i, alpha=-0.6, Fracpoint=1., detect_thresh=5., isl_thresh=3.,  rms_box=(175,43), assocR=None, rootname='model', find_sources_args=""):
        
        self.InMapName = inputmap
        self.InCatName = inputcat
        self.SimI = i
        self.rootname = rootname
        self.find_sources_args = find_sources_args
        
        if not os.path.isdir(dirname):
            os.system("mkdir -p {name:s}".format(name=dirname))
        # initialise names
        self.DirName = dirname
        self.MapName = '{dir:s}/{name:s}{i1:02d}_map.fits'.format(dir=dirname,name=rootname,i1=i)
        self.FlatMapName = '{dir:s}/{name:s}{i1:02d}_flatmap.fits'.format(dir=dirname,name=rootname,i1=i)
        self.InCatName = '{dir:s}/{name:s}{i1:02d}_incat.fits'.format(dir=dirname,name=rootname,i1=i)
        self.OutCatName = '{dir:s}/{name:s}{i1:02d}_outcat.fits'.format(dir=dirname,name=rootname,i1=i)
        self.MatchInOutName = '{dir:s}/{name:s}{i1:02d}_inoutcat.fits'.format(dir=dirname,name=rootname,i1=i)
        self.MatchOutInName = '{dir:s}/{name:s}{i1:02d}_outincat.fits'.format(dir=dirname,name=rootname,i1=i)
        
        # initialise parameters
        self.Nsrc = N
        self.Fmin = Fmin
        self.Fmax = Fmax
        self.FracPoint = Fracpoint
        
        self.alpha = alpha  # flux counts slope
        
        
        # detection 
        self.detect_thresh = detect_thresh
        self.isl_thresh = isl_thresh
        self.rms_box = rms_box
        
        # matching
        # if no assosication radius, get from input image
        if assocR is None:
            infits = pf.open(inputmap)
            inhead = infits[0].header
            beam_maj = inhead.get('BMAJ')
            # or it is not there, check the history
            if beam_maj is None:
                hist = inhead.get_history()
                for i in range(len(hist)):
                    t = hist[i]
                    if 'BMAJ=' in t:
                        tt = t.split()
                        beam_maj = float(tt[tt.index('BMAJ=')+1])
                        break
            fwhm_beam = beam_maj*3600.
            assocR = fwhm_beam*0.25
        self.assocR = assocR
        
        # input files
        if file_exists(residualmap):
            self.ResidMap = residualmap
        else:
            self.ResidMap = None
            
        # input files
        if file_exists(pbmap):
            self.PBMap = pbmap
        else:
            self.PBMap = None
        
        return
    
    def make_model_map(self, clobber=False):
        
        ## for now by default do not clobber ##
        #if not file_exists(self.MapName):
        if (file_exists(self.MapName) and clobber) or (not file_exists(self.MapName)):
            os.system("rm -rf {fle:s}".format(fle=self.MapName))
            
            make_model_map(self.MapName, self.InCatName, self.ResidMap, self.InMapName,  Fmin=self.Fmin, Fmax=self.Fmax, alpha=self.alpha, Nsrc=self.Nsrc, Fracpoint=self.FracPoint)
            
            
        if (file_exists(self.FlatMapName) and clobber) or (not file_exists(self.FlatMapName)):
            os.system("rm -rf {fle:s}".format(fle=self.FlatMapName))
            make_model_flatmap(self.MapName, self.PBMap, self.FlatMapName)
            
            
        return


    def plot_sources_map(self):
        
        
        dirname = self.DirName
        i = self.SimI
        rootname = self.rootname
        inputcat = self.InCatName
        outputcat = self.OutCatName
        
        savename = '{dir:s}/sources_{name:s}{i1:02d}.png'.format(dir=dirname,name=rootname,i1=i)
        
        plot_input_output_cat(inputcat, outputcat, savename)
        
        return


    def find_sources_map(self, cmd, clobber=False):
        
        if (file_exists(self.OutCatName) and clobber) or (not file_exists(self.OutCatName)):
            #find_sources( self.DirName, self.MapName, self.OutCatName, thresh_pix=self.detect_thresh, thresh_isl=self.isl_thresh, rms_box=self.rms_box , clobber=clobber)
            
            os.system('rm -rf {outcat}'.format(outcat=self.OutCatName))
            
            dmap = self.FlatMapName.split('/')[-1]
            mmap = self.MapName.split('/')[-1]
            ddir = self.DirName
            run_cmd = cmd.replace('DETECT', dmap).replace('IMAGE', mmap)
            find_sources_cmd( ddir , run_cmd )
            
            os.system('cp -r {findcat} {outcat}'.format(findcat=self.MapName.replace('.fits','.pybdsm.fits'), outcat=self.OutCatName))
        
        return
    
    def match_catalogues(self, clobber=False):
        
        
        os.system('stilts tcopy in=%s.tempascii out=%s  ifmt=ascii ofmt=fits' %(self.InCatName, self.InCatName))
        
        if (file_exists(self.MatchInOutName) and clobber) or (not file_exists(self.MatchInOutName)):
            os.system('stilts tskymatch2 in1=%s in2=%s out=%s error=%f find=best join=all1' %(self.InCatName, self.OutCatName, self.MatchInOutName, self.assocR))
        
        if (file_exists(self.MatchOutInName) and clobber) or (not file_exists(self.MatchOutInName)):
            os.system('stilts tskymatch2 in1=%s in2=%s out=%s error=%f find=best join=all1' %(self.OutCatName, self.InCatName, self.MatchOutInName, self.assocR))
        
        
        return
    

def rms_stat(RMSMap, fullout=False):
    rms = pf.getdata(RMSMap)
    rmsdat = np.ma.masked_where(np.isnan(rms), rms).compressed()
    rms_mean, rms_sig = clip_mean_sig_nw(rmsdat,3.,5)
    SIG = rms_mean
    SIG_MIN = rmsdat.min()
    SIG_MAX = rmsdat.max()
    SIG_MED = np.median(SIG)

    print 'rms mean : %.3f mJy' %(SIG*1000.)
    if fullout:
        return SIG, SIG_MIN, SIG_MAX, SIG_MED
    else:
        return SIG


def get_differential_area(RMSMap, irange, detect_thresh=5.0):
    rms = pf.getdata(RMSMap)
    rmsdat = np.ma.masked_where(np.isnan(rms), rms).compressed()
    rmshead = pf.getheader(RMSMap)
    # calculate the fraction of area available for sources to be found
    area_i = np.ones(len(irange))
    pixels = np.sum(rmsdat < irange[-1]/detect_thresh)
    ps = rmshead.get('CDELT2') #pixel scale : pix = ps degrees
    degtosterad = (np.pi/180.)**2.
    Area_tot =  pixels*ps*ps * degtosterad
    for i in range( len(irange)-2 ): 
        #incnt = np.sum( (influx>irange[i]) * (influx<irange[i+1]) )
        #cnt = np.sum( (outflux>irange[i]) * (outflux<irange[i+1]) )
        # area of map in which source could be found (S > 5 sigma, sigma < S/5)
        #rmslim = np.sqrt(np.product(irange[i:i+2]))/detect_thresh
        rmslim = irange[i+2]/detect_thresh
        pixels = np.sum(rmsdat < rmslim)
        if pixels == 0:
            area_i[i] = np.nan
            continue
        Area = pixels*ps*ps
        Area_st = Area * degtosterad
        area_i[i] = Area_st/Area_tot
    return area_i

def Completeness(MatchedInOutCatList, SIG, irange, area_i, detect_thresh=5., showall=False, savename='', title='', plotvline=None, rmsfits=None):
    ####### COMPLETENESS ################
    print 'CALCULATE COMPLETENESS'



    #fraction_S = []
    #snr_S = []
    #inpeak_all = []
    #outpeak_all = []
    #outflux_all = []
    #rms_all=[]
    fraction_total_MC = []
    fraction_total_MC_n = []
    fraction_total_MC_d = []

    # count sources
    for MatchedInOutCat in MatchedInOutCatList:
        #if verbose : print 'RUN %i' %(i)
        print MatchedInOutCat
        matched_cat_list_inout = load_fits(MatchedInOutCat)

        
        inpeak = matched_cat_list_inout.Peak_flux_1.copy()
        outpeak = matched_cat_list_inout.Peak_flux_2.copy()
        influx = matched_cat_list_inout.Total_flux_1.copy()
        outflux = matched_cat_list_inout.Total_flux_2.copy()
        rms = matched_cat_list_inout.Isl_rms.copy()
        detect = np.isfinite(matched_cat_list_inout.Isl_rms)
        nsrc = len(matched_cat_list_inout)
        
        
        mask_high = outpeak/inpeak > 1.3
        mask_low = outpeak/inpeak < 0.77
        mask_bad_m = np.ma.mask_or(mask_low, mask_high)
        mask_bad = np.where(mask_bad_m)[0] 
        outflux[mask_bad] = np.nan*np.ones(len(mask_bad))
        outpeak[mask_bad] = np.nan*np.ones(len(mask_bad))
        
        total_dfrac = np.ones(len(irange)-1)*np.nan
        total_dfrac_n = np.ones(len(irange)-1)*np.nan
        total_dfrac_d = np.ones(len(irange)-1)*np.nan
        for k in range(len(irange)-1):
            mask = 1-((influx > irange[k]) * (influx < irange[k+1]))
            influx_mask = np.ma.masked_where(mask, influx).compressed()
            nsrc_mask = len(influx_mask)
            outflux_mask = np.ma.masked_where(mask, outflux).compressed()
            detections = np.sum(outflux_mask > 0.)  #-1 are non-detections
            #detections = np.sum(np.ma.masked_where(outpeak/rms >pNrange[k],detect).compressed())
            #print prange[k],nsrc_mask
            if nsrc_mask == 0.:
                total_dfrac[k] = 0.
            else:
                total_dfrac[k] = float(detections)/nsrc_mask
                total_dfrac_n[k] = float(detections)
                total_dfrac_d[k] = nsrc_mask
                #if verbose: 
                    #print  '%.2f %.2f %.2f %.2f' %(irange[k]*1000,  irange[k+1]*1000, float(detections), nsrc_mask)
            
        fraction_total_MC.append(total_dfrac)
        fraction_total_MC_n.append(total_dfrac_n)
        fraction_total_MC_d.append(total_dfrac_d)
        outpeakd = np.ma.masked_where(detect==0,np.array(outpeak)).compressed()
        outfluxd = np.ma.masked_where(detect==0,np.array(outflux)).compressed()
        inpeakd = np.ma.masked_where(detect==0,np.array(inpeak)).compressed()
        influxd = np.ma.masked_where(detect==0,np.array(influx)).compressed()
        rmsd = np.ma.masked_where(detect==0, np.array(rms)).compressed()
        detectd =  np.ma.masked_where(detect==0,detect).compressed()

        #inpeak_all.append(inpeak)
        #outpeak_all.append(outpeak)
        #rms_all.append(rms)
        #outflux_all.append(outflux)
    
        print 'FLUX        DF Tot Frac'
        for k in range(len(total_dfrac)):                                     
            print '%10.3f %3.0f %3.0f %5.3f' %( irange[k]*1000., total_dfrac_n[k], total_dfrac_d[k], total_dfrac_n[k]/total_dfrac_d[k] )

    irangep = (irange[:-1] + irange[1:])/2.

    
    
    #inpeak_all = np.array(inpeak_all)
    #outpeak_all = np.array(outpeak_all)
    #outflux_all = np.array(outflux_all)
    #rms_all = np.array(rms_all)
    
    #SNRl = outpeak_all/rms_all
    #Sratl = outflux_all/outpeak_all

    #SNR = SNRl.flatten() #[0]
    #Srat = Sratl.flatten() #[0]

    #sys.exit()
        

    
    import matplotlib.pyplot as plt

    #if mode == 'e':
    if 1:
        fraction_total_MC_n = np.array(fraction_total_MC_n)
        fraction_total_MC_d = np.array(fraction_total_MC_d)

        fraction_total_MC = fraction_total_MC_n/fraction_total_MC_d
        #fraction_total_MC_mean = np.average(fraction_total_MC, axis=0)
        #fraction_total_MC_sig = np.std(fraction_total_MC, axis=0)

        fraction_total_MC_n_av = np.average(fraction_total_MC_n, axis=0)
        fraction_total_MC_n_sig = np.std(fraction_total_MC_n, axis=0)
        fraction_total_MC_d_av = np.average(fraction_total_MC_d, axis=0)
        fraction_total_MC_d_sig = np.std(fraction_total_MC_d, axis=0)

        #fraction_total_MC_mean = fraction_total_MC_n_av/fraction_total_MC_d_av
        #fraction_total_MC_sig = fraction_total_MC_mean*np.sqrt( (fraction_total_MC_n_sig/fraction_total_MC_n_av)**2. + (fraction_total_MC_d_sig/fraction_total_MC_d_av)**2. )
        fraction_total_MC_mean = np.average(fraction_total_MC_n/fraction_total_MC_d, axis=0)
        fraction_total_MC_sig = np.std(fraction_total_MC_n/fraction_total_MC_d, axis=0)
        
        
        if rmsfits is not None:
        
            print 'getting areas'
            #rmsdat = 0.79*pf.getdata(rmsfits)
            rmsdat = pf.getdata(rmsfits)
            rmsdat = rmsdat[np.isfinite(rmsdat)]
            rmsdat = rmsdat[rmsdat>0]
            rmsdat = rmsdat.flatten()
            #rmsdat = rmsdat/corFlux
            #rmshead = pf.getheader(rmsfits)
            #ps = rmshead.get('CDELT2')

            Area_rms_int = np.nan*np.ones(len(irangep))
            pixel_tot = 1.0*np.sum(np.isfinite(rmsdat))
            for i in range(len(irangep)):
                fp = irange[i+1]
                rmslim = fp/5.
                pixel = np.sum(rmsdat < rmslim)
                #area = pixel *ps*ps
                #Area_rms_int[i] = area/Atot
                Area_rms_int[i] = pixel/pixel_tot
            #Area_rms_int[pixels==0] = np.nan
        
        
        #f1 = plt.figure(figsize = (7.0,4.5))
        #ax1 = plt.subplot(111)
        f1, ax1 = paper_single_ax()
        ax1.minorticks_on()
        #for i in range(len(fraction_total_MC)):
        #ax1.semilogx(irange,fraction_total_MC[i])
        
        
        ax1.semilogx(1000.*irangep,fraction_total_MC_mean,'k', label='Fraction')
        ax1.semilogx(1000.*irangep,fraction_total_MC_mean+fraction_total_MC_sig,'k:')
        ax1.semilogx(1000.*irangep,fraction_total_MC_mean-fraction_total_MC_sig,'k:')
        if rmsfits is not None:
            #ax1.semilogx(1000.*irangep,Area_rms_int,'b')
            
            ax1.semilogx(1000.*irangep,fraction_total_MC_mean/Area_rms_int,'r', label='Area-corrected')
            ax1.semilogx(1000.*irangep,fraction_total_MC_mean/Area_rms_int+fraction_total_MC_sig,'r:')
            ax1.semilogx(1000.*irangep,fraction_total_MC_mean/Area_rms_int-fraction_total_MC_sig,'r:')
            
        ax1.vlines(3.*SIG*1000.,0,1.05,linestyle='dashed')
        if plotvline is not None:
            ax1.vlines(plotvline*1000,0,1.05,linestyle='dotted')
        ax1.set_ylabel('Detected Fraction')
        ax1.set_xlabel('$S_{\mathrm{int}}$ [mJy] ')
        ax1.set_ylim(0,1.05)
        ax1.set_xlim(1000.*irangep.min(), 1000.*irangep.max())
        #ax1.set_xlim(5., 100.)
        ax1.minorticks_on()
        plt.subplots_adjust(bottom=0.15,right=0.95)
        ax1.set_title(title)
        ax1.legend(loc='lower right')
        if showall:
            ax1.semilogx(1000.*irangep,fraction_total_MC.transpose())
        #f1.savefig('%s_Detection_fraction_int_i.png' %(title))
        #f1.savefig('%s_Detection_fraction_int_i.eps' %(title),dpi=600)
        #ax1.yaxis.set_major_locator(MultipleLocator(1.0))
        #ax1.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
        ax1.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
        ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x) if x < 1 else ('%.0f')%(x)))
        #ax1.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
        fig_save_many(f1, '%s_Detection_fraction_int_i' %(savename),dpi=600)
        np.save('%s_Detection_fraction_int_i.npy' %(savename), [area_i,1000.*irangep,fraction_total_MC_mean,fraction_total_MC_sig,fraction_total_MC_mean/Area_rms_int])

    #if mode == 'e':
    if 1:
        #f1 = plt.figure(figsize = (7.0,4.5))
        #ax1 = plt.subplot(111)
        f1, ax1 = paper_single_ax()
        ax1.minorticks_on()

        #Completeness is sum of sources detected above flux S over sum of sources which we could have detected (ie input distribution scaled down by area fraction)
        completeness_diff_n = (fraction_total_MC_n) #/(area_i)
        completeness_diff_d = (fraction_total_MC_d)  #*(area_i)
        
        print completeness_diff_d.shape
        print completeness_diff_n.shape
        
        nanind = np.where(np.isnan(completeness_diff_d) & np.isnan(completeness_diff_n))
        dennanind = np.where(np.isnan(completeness_diff_d) & np.isfinite(completeness_diff_n))
        numnanind = np.where(np.isnan(completeness_diff_n) & np.isfinite(completeness_diff_n))
        
        completeness_diff_d[nanind] = 1
        completeness_diff_n[nanind] = 1
        
        #completeness_diff_d = completeness_diff_d[np.isfinite(completeness_diff_n)]
        #completeness_diff_n = completeness_diff_n[np.isfinite(completeness_diff_n)]
        ##irangep = irangep[np.isfinite(completeness_diff_n)]
        #completeness_diff_n = completeness_diff_n[np.isfinite(completeness_diff_d)]
        #completeness_diff_d = completeness_diff_d[np.isfinite(completeness_diff_d)]
        #irangep = irangep[np.isfinite(completeness_diff_d)]

        print completeness_diff_d.shape
        print completeness_diff_n.shape

        completeness_int = np.zeros(completeness_diff_n.shape)
        if rmsfits is not None:
            completeness_int_area = np.zeros(completeness_diff_n.shape)
        for i in range(len(irange)-1):
            completeness_int[:,i] = np.sum(completeness_diff_n[:,i:],axis=1)/np.sum(completeness_diff_d[:,i:],axis=1)
            #if rmsfits is not None:
                #completeness_int_area[:,i] = np.sum(completeness_diff_n[:,i:]/Area_rms_int,axis=1)/np.sum(completeness_diff_d[:,i:]/Area_rms_int,axis=1)
        
        Comp = np.nanmean(completeness_int, axis=0)
        Comp_sd = np.nanstd(completeness_int, axis=0)
        #if rmsfits is not None:
            #Comp_Area = np.nanmean(completeness_int_area, axis=0)
            #Comp_Area_sd = np.nanstd(completeness_int_area, axis=0)

        #if rmsfits is not None:
            #ax1.semilogx(1000.*irangep,Comp_Area,'r')
            #ax1.semilogx(1000.*irangep,Comp_Area+Comp_Area_sd,'r:')
            #ax1.semilogx(1000.*irangep,Comp_Area-Comp_Area_sd,'r:')
        try:
            ax1.semilogx(1000.*irangep,Comp,'k')
        except:
            print irangep
            print Comp
            return
        
        
        ax1.semilogx(1000.*irangep,Comp+Comp_sd,'k:')
        ax1.semilogx(1000.*irangep,Comp-Comp_sd,'k:')

        ax1.vlines(3.*SIG*1000.,0,1.05,linestyle='dashed')
        if plotvline is not None:
            print plotvline
            ax1.vlines(plotvline*1000,0,1.05,linestyle='dotted')
        ax1.set_ylabel('Completeness ($>S_{\mathrm{int}}$)')
        ax1.set_xlabel('$S_{\mathrm{int}}$ [mJy] ')
        ax1.set_ylim(0,1.05)
        ax1.set_xlim(1000.*irangep.min(), 1000.*irangep.max())
        #ax1.set_xlim(5., 100.)
        ax1.minorticks_on()
        plt.subplots_adjust(bottom=0.15,right=0.95)
        ax1.set_title(title)
        #f1.savefig('%s_Completeness_int_i.png' %(title))
        #f1.savefig('%s_Completeness_int_i.eps' %(title),dpi=600)
        #ax1.yaxis.set_major_locator(MultipleLocator(1.0))
        #ax1.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
        ax1.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
        ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x) if x < 1 else ('%.0f')%(x)))
        #ax1.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
        fig_save_many(f1, '%s_Completeness_int_i' %(savename),dpi=600 )
        np.save('%s_Completeness_int_i.npy' %(savename), [1000.*irangep,Comp])

        for i in range(len(Comp)):
            if np.isnan(Comp[i]): Comp[i] = 0.
        #x = np.logspace(i1,i2,Npoints*10)
        x = np.logspace(irangep.min(),irangep.max(),len(irangep)*10)
        
        y = np.interp(x, irangep, Comp)
        print '#Completeness#'
        print '90 per cent complete %.2f mJy' %(x[np.sum(y<0.9) -1]*1000)
        print '95 per cent complete %.2f mJy' %(x[np.sum(y<0.95) -1]*1000)
        print '97.5 per cent complete %.2f mJy' %(x[np.sum(y<0.975) -1]*1000)


def invert_map(input_map_name, invert_map_name):
    
    print "inverting map"
    
    infits = pf.open(input_map_name)
    inhead = infits[0].header
    input_data = infits[0].data
    inverted_data = -1.*input_data
    
    pf.writeto(invert_map_name, inverted_data, inhead, clobber=True)
    return

def Reliability(PBmap, RESIDUALmap, SIG, invert_map_dir, invert_map_name, SOURCElist, irange, findcmd, detect_thresh=5.0, isl_thresh=3.0, clobber_find=True, title='', savename='', rms_box=(175,43), plotvline=None ):
    #FIELD = 'A2256'
    #target = FIELD
    #RESIDUALmap = '/data1/wwilliams/phd/gmrt_bootes/a2256_mosaic/fits/%s.MOSAIC.MASKED.RES.FITS' %(FIELD)
    #A2256.MOSAIC.3RES.RES.FITS
    #invert_map_dir = '/data1/wwilliams/phd/gmrt_bootes/a2256_mosaic/fits/invert_maps/'
    #SOURCElist = '/data1/wwilliams/phd/gmrt_bootes/a2256_mosaic/fits/%s.MOSAIC.MASKED.pybdsm.clean.fits' %(FIELD)
    #TITLE = FIELD
#if 1:
    print 'CALCULATE RELIABILITY'


    if not os.path.isdir(invert_map_dir):
        os.system("mkdir -p {name:s}".format(name=invert_map_dir))

    
    irangep = (irange[:-1] + irange[1:])/2.

    #t=pf.open(SOURCElist)[1].data
    #NSOURCES = t.shape[0]
    #FLUXMAX = t.Total_flux.max()
    #FLUXMIN = t.Total_flux.min()

    #detect_thresh = 5.0
    #isl_thresh = 3.0
    #clobber_find = 0

    # detect fake sources
    #invert_map_name = '%s_RES_INVERTED.FITS' %(title)

    flat_invert_map_name = invert_map_name.replace('fits', 'flat.fits')
    if not os.path.exists(invert_map_dir+invert_map_name):
        invert_map(RESIDUALmap, invert_map_dir+invert_map_name)
    if not os.path.exists(invert_map_dir+flat_invert_map_name):
        make_model_flatmap(invert_map_dir+invert_map_name, PBmap, invert_map_dir+flat_invert_map_name) 
    #invert_map(INPUTmap, invert_map_dir+invert_map_name)

    outcat = invert_map_dir+invert_map_name.replace('.fits','.pybdsm.fits')
    
    exists = os.path.isfile(outcat)
    if (exists and clobber_find) or (not exists):        
        #find_sources( invert_map_dir, invert_map_name, outcat, thresh_pix=detect_thresh, thresh_isl=isl_thresh, rms_box=rms_box )
        
        dmap = flat_invert_map_name
        mmap = invert_map_name
        run_cmd = findcmd.replace('DETECT', dmap).replace('IMAGE', mmap)
        find_sources_cmd(invert_map_dir, run_cmd)
        
    invert_cat = load_fits(outcat)
    real_cat = load_fits(SOURCElist)

    # remove some artefact detections
    #bad_src_list = [12, 13, 14]
    #bad_src_list = [0, 31, 11, 41, 15, 28]
    bad_src_list = []
    good_ind = [ i for i,src in enumerate(invert_cat['Source_id']) if src not in bad_src_list ]
    invert_cat = invert_cat[good_ind]

    #false_fluxes = invert_cat['Total_flux']
    #real_fluxes = real_cat['Total_flux']
    false_fluxes = invert_cat['Peak_flux']
    real_fluxes = real_cat['Peak_flux']


    def countinbins(xbins, ydata):
        ycounts = np.array([ np.sum((ydata > xbins[xi]) & (ydata <= xbins[xi+1]) )  for xi in range(len(xbins) -1) ], dtype=float )
        return ycounts


    def countabovebins(xbins, ydata):
        ycounts = np.array([ np.sum(ydata > xbins[xi])   for xi in range(len(xbins) -1) ], dtype=float )
        return ycounts

    import matplotlib.pyplot as plt

    FD_n = countinbins(irange, false_fluxes)
    RD_n = countinbins(irange, real_fluxes)
    RD_n_mask = RD_n.copy()
    FD_n_mask = FD_n.copy()
    FD_n_mask[RD_n_mask==0] = np.nan
    RD_n_mask[RD_n_mask==0] = np.nan
    RD_n_mask[FD_n_mask==0] = np.nan
    FD_n_mask[FD_n_mask==0] = np.nan
    FDR_n = FD_n/RD_n_mask
    print FD_n
    print RD_n
    print FD_n_mask
    print RD_n_mask
    #sFDR_n = FDR_n*np.sqrt(1./FD_n_mask + 1./RD_n_mask)
    sFDR_n = np.sqrt(FD_n_mask/RD_n_mask**2. + FD_n_mask**2./RD_n_mask**3.)
    print sFDR_n

    #f1 = plt.figure(figsize = (7.0,4.5))
    #ax1 = plt.subplot(111)
    f1, ax1 = paper_single_ax()
    ax1.minorticks_on()
    ax1.set_xscale('log')
    #for i in range(len(Rfraction_MC)):
    #ax1.semilogx(irange,Rfraction_MC[i])
    #print i, average(Rfraction_MC[i])
    ax1.semilogx(1000.*irangep,FDR_n,'k')
    ax1.semilogx(1000.*irangep,FDR_n+sFDR_n,'k:')
    ax1.semilogx(1000.*irangep,FDR_n-sFDR_n,'k:')
    #print '90 per cent reliable above integrated flux of %.4f' %(irange[np.sum(abs(Rfraction_MC_mean)>0.10)-1]*1000.)
    #print '95 per cent reliable above integrated flux of %.4f' %(irange[np.sum(abs(Rfraction_MC_mean)>0.05)-1]*1000.)
    ax1.vlines(3.*SIG*1000.,0,1.05,linestyle='dashed')
    if plotvline is not None:
        ax1.vlines(plotvline*1000,0,1.05,linestyle='dotted')
    ax1.set_ylabel('False Detection Rate')
    ax1.set_xlabel('$S_{\mathrm{int}}$ [mJy] ')
    ax1.set_ylim(0,1.05)
    ax1.set_xlim(1000.*irangep.min(), 1000.*irangep.max())
    #ax1.set_xlim(5., 100.)
    ax1.minorticks_on()
    plt.subplots_adjust(bottom=0.15,right=0.95)
    ax1.set_title(title)
    #f1.savefig('%s_Rely_fraction_1sig_i.png' %(title))
    #f1.savefig('%s_Rely_fraction_1sig_i.eps' %(title),dpi=600)
    #ax1.yaxis.set_major_locator(MultipleLocator(1.0))
    #ax1.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    ax1.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x) if x < 1 else ('%.0f')%(x)))
    #ax1.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    fig_save_many(f1, '%s_Rely_fraction_1sig_i' %(savename),dpi=600)
    ax1.set_ylim(0,0.41)
    fig_save_many(f1, '%s_Rely_fraction_1sig_i_zoom' %(savename),dpi=600)

    np.save('%s_Rely_fraction_1sig_i.npy' %(savename), [1000.*irangep,FDR_n,sFDR_n])

    FD = [np.sum(FD_n[i:]) for i in range(len(FD_n)) ]
    RD = [np.sum(RD_n[i:]) for i in range(len(RD_n)) ]
    sFD = [np.sqrt(np.sum(np.sqrt(FD_n[i:]))) for i in range(len(FD_n))  ]
    sRD = [np.sqrt(np.sum(np.sqrt(RD_n[i:]))) for i in range(len(RD_n))  ]
    FD = countabovebins(irange, false_fluxes)
    RD = countabovebins(irange, real_fluxes)
    RD_mask = RD.copy()
    RD_mask[RD_n==0] = np.nan
    FD_mask = FD.copy()
    FD_mask[FD==0] = np.nan
    #sFD = np.sqrt(np.sum([np.sqrt(FD_ni) for FD_ni in FD_mask]))
    #sRD = np.sqrt(np.sum([np.sqrt(RD_ni) for RD_ni in RD_mask]))
    Rely = 1. - FD/RD_mask
    #Rely_sd = Rely*np.sqrt(sFD/FD_mask + sRD/RD_mask)
    Rely_sd = np.sqrt( (FD/(RD**2.)) + ( (FD**2.)/(RD**3.)))
    #Rely_sd = Rely*np.sqrt(1./(RD_mask-FD_mask) + 1./RD_mask)

    #f1 = plt.figure(figsize = (7.0,4.5))
    #ax1 = plt.subplot(111)
    f1,ax1 = paper_single_ax()
    ax1.set_xscale('log')
    ax1.minorticks_on()
    ax1.plot(1000.*irangep,Rely,'k')
    ax1.semilogx(1000.*irangep,Rely,'k')
    ax1.semilogx(1000.*irangep,Rely+Rely_sd,'k:')
    ax1.semilogx(1000.*irangep,Rely-Rely_sd,'k:')
    ax1.vlines(3.*SIG*1000.,0,1.05,linestyle='dashed')
    if plotvline is not None:
        ax1.vlines(plotvline*1000,0,1.05,linestyle='dotted')
    ax1.set_ylabel('Reliability ($>S_{\mathrm{int}}$)')
    ax1.set_xlabel('$S_{\mathrm{int}}$ [mJy] ')
    ax1.set_ylim(0,1.05)
    ax1.set_xlim(1000.*irangep.min(), 1000.*irangep.max())
    #ax1.set_xlim(5., 100.)
    ax1.set_ylim(0., 1.05)
    ax1.minorticks_on()
    plt.subplots_adjust(bottom=0.15,right=0.95)
    ax1.set_title(title)
    #f1.savefig('%s_Rely_1sig_i.png' %(title))
    #f1.savefig('%s_Rely_1sig_i.eps' %(title),dpi=600)
    #ax1.yaxis.set_major_locator(MultipleLocator(1.0))
    #ax1.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    ax1.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x) if x < 1 else ('%.0f')%(x)))
    #ax1.yaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
    fig_save_many(f1, '%s_Rely_1sig_i' %(savename),dpi=600)
    ax1.set_ylim(0.55, 1.04)
    fig_save_many(f1, '%s_Rely_1sig_i_zoom' %(savename),dpi=600)
    np.save('%s_Rely_1sig_i.npy' %(savename), [1000.*irange,Rely])


    print 'FLUX        Rely FD RD'
    for k in range(len(Rely)):                                     
        print '%10.3f %5.3f %3.0f %3.0f %4.2f' %( irange[k]*1000., Rely[k], FD_n[k], RD_n[k], FD_n[k]/RD_n[k] )


    for i in range(len(Rely)):
        if np.isnan(Rely[i]): Rely[i] = 0.
    x = np.logspace(irange.min(),irange.max(),len(irange)*10)
    y = np.interp(x, irangep, Rely)
    print '#Reliability#'
    print '90 per cent reliability above flux of %.2f mJy' %(x[np.sum(y<0.90) -1]*1000)
    print '95 per cent reliability above flux of %.2f mJy' %(x[np.sum(y<0.95) -1]*1000)
    print '97.5 per cent reliability above flux of %.2f mJy' %(x[np.sum(y<0.975) -1]*1000)

    return




#def main():
    #return

