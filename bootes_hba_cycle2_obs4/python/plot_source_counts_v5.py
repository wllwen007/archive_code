
import matplotlib as mpl
mpl.use('Agg')
import pylab as pl
import numpy as np
from fits_util import *
import scipy
from scipy.optimize import leastsq
import warnings, os
warnings.filterwarnings('ignore')
import plot_util as pp
import argparse
from sky_util import angular_separation
import astropy.coordinates as ac


#import matplotlib as mpl
import matplotlib as mpl
mpl.rc_file('/home/wwilliams/.config/matplotlib/matplotlibrc')  # <-- the file containing your settings

#mpl.rcParams["text.usetex"] = False
#mpl.rcParams["font.family"] = "sans-serif"


# v3
# add multiple source correction

def time_smearing2(delta_T,resolution,delta_Theta):
        #Same as above but provides the flux loss for a given time averaging
        Reduction = 1-1.22E-9*(delta_Theta/resolution)**2.0 * delta_T**2.0
        #print 'At radius %s the source will have %s percent of its flux if data smoothed to %s'%(delta_Theta,Reduction,delta_T)

        return Reduction


def bandwidth_smearing2(delta_freq,freq,resolution,delta_Theta):
        # Same as above but gives the flux loss for a given frequency averaging.

    beta = (delta_freq/freq) * (delta_Theta/resolution)
    gamma = 2*(np.log(2)**0.5)
    Reduction = ((np.pi**0.5)/(gamma * beta)) * (scipy.special.erf(beta*gamma/2.0))

    #print 'At radius %s a source will only have %s percent of its flux if data smoothed in freq to %s'%(delta_Theta,Reduction,delta_freq)

    return Reduction


def Tpoly(c,x):
    """ y = Sum { c(i)*x^i }, i=0,len(c)"""
    y = np.zeros(len(x))
    for i in range(len(c)):
        y += c[i]*(x**i)
    return y

  
def  get_IC10_counts():
  fluxdat = pl.loadtxt('/local/wwilliams/phd/gmrt_bootes/bootes_analysis/IC10_counts.dat').transpose()
  
  bin_centres = fluxdat[0]/1000.
  counts = fluxdat[1]
  count_errors = (fluxdat[4]+fluxdat[5])/2.
  frq = 150.
  return frq, bin_centres ,counts, count_errors
  
def  get_skads_counts(ttype=[[0,1,2],[0,1,2,3,4]], freq=150., f=150.):
  print 'SKADS'
  sft = ['None', 'Quiescent', 'Starburst']
  agnt = ['None', 'RQ',' FRI','FRII','GPS']
  
  sftype = ttype[0]
  agntype = ttype[1]
  
  selsft = [ sft[i] for i in sftype ]
  selagnt = [ agnt[i] for i in agntype ]
  print ' agn : %s : %s ' %(agntype, ','.join(selagnt))  #None [0] Radio-Quiet [1] FRI [2] FRII [3] GPS [4]
  print ' sf : %s : %s' %(sftype, ','.join(selsft))   #None [0] Quiescent [1] Starburst [2]
  
  data = load_fits('/local/wwilliams/projects/skads/wilman_cat_all.fits')
  
  z = data['redshift']
  sftypes =data['sftype']
  agntypes = data['agntype']
  if freq == 150.:
    ff = data['itot_151']
  elif freq == 610.:
    ff = data['itot_610']
  
  
  
  m = pl.zeros(len(agntypes),dtype=int)
  for agnt in agntype:
      m += (agntypes == agnt)
  for sft in sftype:
      m += (sftypes == sft)
  mask = m>0
  
  #mask = pl.zeros(len(agntypes),dtype=bool)
  #for i in range(len(agntypes)):
    #if agntypes[i] in agntype:
      #mask[i] = 1
    #if sftypes[i] in sftype:
      #mask[i] = 1
  fluxes = 10.**ff
  flux = fluxes[mask]
  #flux = pl.ma.masked_where(mask, fluxes).compressed()
  print ' n sources = %i' %(len(flux))
  #b1 = 1.e-4
  #b2 = 4.
  #bins = pl.logspace(pl.log10(b1),pl.log10(b2),35.)
  #bin_centres = pl.zeros( len(bins) - 1 )
  #for i in range(len(bins) - 1):
    #bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  if freq == 150.:
    b1 = 1.e-4
    b2 = 4
  elif freq == 610.:
    b1 = 0.5e-4
    b2 = 6

  A = 51.18
  beta = -1.436
  Nperbin = 100.
  bins=[b1]
  while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
  bin_centres = pl.zeros( len(bins) - 1 )
  for i in range(len(bins) - 1):
    bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  nbins =  len(bin_centres) 
  bincount = pl.zeros(nbins)
  #dnbins = pl.zeros( nbins )
  counts = pl.zeros( nbins )
  count_errors = pl.zeros( nbins )
  #print bins
  #print nbins
  #print len(bins)
  #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('BinL', 'BinC', 'BinW', 'Count', 'Count/3', 'Norm Count', 'Area deg^2', 'Area st') 
  Area = 20.*20.
  degtosterad = (pl.pi/180.)**2.
  Area_st = Area * degtosterad
  for i in range( nbins ): 
    cnt = pl.sum( (flux>bins[i]) * (flux<bins[i+1]) )
    # area of sim
    bincount[ i ] = cnt
    nbin_err = pl.sqrt( cnt )
    dnbin = bins[i+1] - bins[i]
    counts[ i ] = ((bincount[i]/dnbin)/Area_st)*(bin_centres[i]**2.5)
    count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
    #print '%.4f\t%.4f\t%.4f\t%i\t%i\t%.3f\t%.2f\t%.5f' %(bins[i], bin_centres[ i ], dnbin, cnt, cnt/3., counts[i], Area, Area_st) 
  
  
  
  if freq == 610:
    # scale to required freq
    bin_centres = bin_centres*(f/freq)**alpha
    counts = counts*(f/freq)**(1.5*alpha)
    count_errors = count_errors*(f/freq)**(1.5*alpha)
  
  
  return bin_centres ,counts, count_errors


def  get_skads_counts_sfonly(majorlim=0, less=True, freq=150., f=150.):
  print 'SKADS SF'
    
  data = load_fits('/local/wwilliams/projects/skads/wilman_cat_sf_comp.fits')
  
  z = data['redshift']
  #ff = data['itot_151']
  if freq == 150.:
    ff = data['itot_151']
  elif freq == 610.:
    ff = data['itot_610']
  
  major = data['major_axis']
  
  if less:
    mask = np.where(major<majorlim)
  else:
    mask = np.where(major>=majorlim)
      
  fluxes = 10.**ff
  flux = fluxes[mask]
  
  print ' n sources = %i' %(len(flux))
  #b1 = 1.e-4
  #b2 = 4.
  #bins = pl.logspace(pl.log10(b1),pl.log10(b2),35.)
  #bin_centres = pl.zeros( len(bins) - 1 )
  #for i in range(len(bins) - 1):
    #bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  if freq == 150.:
    b1 = 1.e-4
    b2 = 4
  elif freq == 610.:
    b1 = 0.5e-4
    b2 = 6

  A = 51.18
  beta = -1.436
  Nperbin = 200.
  bins=[b1]
  while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
  bin_centres = pl.zeros( len(bins) - 1 )
  for i in range(len(bins) - 1):
    bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  nbins =  len(bin_centres) 
  bincount = pl.zeros(nbins)
  #dnbins = pl.zeros( nbins )
  counts = pl.zeros( nbins )
  count_errors = pl.zeros( nbins )
  #print bins
  #print nbins
  #print len(bins)
  #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('BinL', 'BinC', 'BinW', 'Count', 'Count/3', 'Norm Count', 'Area deg^2', 'Area st') 
  Area = 20.*20.
  degtosterad = (pl.pi/180.)**2.
  Area_st = Area * degtosterad
  for i in range( nbins ): 
    cnt = pl.sum( (flux>bins[i]) * (flux<bins[i+1]) )
    # area of sim
    bincount[ i ] = cnt
    nbin_err = pl.sqrt( cnt )
    dnbin = bins[i+1] - bins[i]
    counts[ i ] = ((bincount[i]/dnbin)/Area_st)*(bin_centres[i]**2.5)
    count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
    #print '%.4f\t%.4f\t%.4f\t%i\t%i\t%.3f\t%.2f\t%.5f' %(bins[i], bin_centres[ i ], dnbin, cnt, cnt/3., counts[i], Area, Area_st) 
  
  
  if freq == 610:
    # scale to required freq
    bin_centres = bin_centres*(f/freq)**alpha
    counts = counts*(f/freq)**(1.5*alpha)
    count_errors = count_errors*(f/freq)**(1.5*alpha)
  
  return bin_centres ,counts, count_errors

def  get_skads_counts_agnonly(freq=150., f=150.):
  print 'SKADS AGN'
    
  data = load_fits('/local/wwilliams/projects/skads/wilman_cat_agn.fits')
  
  z = data['redshift']
  if freq == 150.:
    ff = data['itot_151']
  elif freq == 610.:
    ff = data['itot_610']
  
  
  #mask = np.where(major>majorlim)
  flux = 10.**ff
  #flux = fluxes[mask]
  
  print ' n sources = %i' %(len(flux))
  #b1 = 1.e-4
  #b2 = 4.
  #bins = pl.logspace(pl.log10(b1),pl.log10(b2),35.)
  #bin_centres = pl.zeros( len(bins) - 1 )
  #for i in range(len(bins) - 1):
    #bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  if freq == 150.:
    b1 = 1.e-4
    b2 = 4
  elif freq == 610.:
    b1 = 0.5e-4
    b2 = 6
      

  A = 51.18
  beta = -1.436
  Nperbin = 100.
  bins=[b1]
  while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
  bin_centres = pl.zeros( len(bins) - 1 )
  for i in range(len(bins) - 1):
    bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  nbins =  len(bin_centres) 
  bincount = pl.zeros(nbins)
  #dnbins = pl.zeros( nbins )
  counts = pl.zeros( nbins )
  count_errors = pl.zeros( nbins )
  #print bins
  #print nbins
  #print len(bins)
  #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('BinL', 'BinC', 'BinW', 'Count', 'Count/3', 'Norm Count', 'Area deg^2', 'Area st') 
  Area = 20.*20.
  degtosterad = (pl.pi/180.)**2.
  Area_st = Area * degtosterad
  for i in range( nbins ): 
    cnt = pl.sum( (flux>bins[i]) * (flux<bins[i+1]) )
    # area of sim
    bincount[ i ] = cnt
    nbin_err = pl.sqrt( cnt )
    dnbin = bins[i+1] - bins[i]
    counts[ i ] = ((bincount[i]/dnbin)/Area_st)*(bin_centres[i]**2.5)
    count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
    #print '%.4f\t%.4f\t%.4f\t%i\t%i\t%.3f\t%.2f\t%.5f' %(bins[i], bin_centres[ i ], dnbin, cnt, cnt/3., counts[i], Area, Area_st) 
  
  if freq == 610:
    # scale to required freq
    bin_centres = bin_centres*(f/freq)**alpha
    counts = counts*(f/freq)**(1.5*alpha)
    count_errors = count_errors*(f/freq)**(1.5*alpha)
  
  return bin_centres ,counts, count_errors



def get_deZotti_counts(f, alpha=-0.8, value='all'):
    dat = np.genfromtxt('/local/wwilliams/projects/src_counts/deZotti2010_compilation.dat', skip_header=22, dtype=[('logS','<f8'),('cnt','<f8'),('u_e_cnt','<f8'),('l_e_cnt','<f8'),('var','<f8'),('ref','S5')])
    # col. 1 =      log(S [Jy])
# col. 2 =      S^2.5dN/dS[Jy^1.5/sr]
# col. 3-4 =    positive and negative errors on the counts.
# col. 5 =      For surveys covering less than 25 deg^2 the contribution to such errors
#               due to the sampling variance [eq.(6) in De Zotti et al. (2009)];
#               such contributions are negligible for larger areas.
# col. 6 =      references (see the list below)
#    bo08: 2008ApJ...681.1129B
#    br72: 1972AJ.....77..405B
#    ci99: 1999MNRAS.302..222C
#    fo06: 2006ApJS..167..103F
#    gr99: 1999MNRAS.305..297G
#    ho03: 2003AJ....125..465H
#    ke08: 2008ApJS..179...71K
#    ow08: 2008AJ....136.1889O
#    ri00: 2000ApJ...533..611R
#    se08A:2008MNRAS.386.1695S
#    wi97: 1997ApJ...475..479W
    
    if value != 'all':
        dat = dat[dat['ref']==value]
    
    bin_centres = 10**(dat['logS'])  # to Jy
    counts = dat['cnt']
    count_errors = np.array([dat['l_e_cnt'], dat['u_e_cnt']])
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors

def get_Ciliegi_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/ciliegi_1999.dat').transpose()
    bin_centres = dat[2]/1e3  # to Jy
    counts = dat[5]
    count_errors = dat[6]
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors

def get_Richards_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/richards_2000.dat').transpose()
    bin_centres = dat[2]/1e6  # to Jy
    counts = dat[6]
    count_errors = dat[7]
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors


def get_Seymour_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/seymour_2004.dat').transpose()
    bin_centres = dat[2]/1e6  # to Jy
    counts = dat[5]
    count_errors = dat[6]
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors

def get_Prandoni_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/prandoni_2001.dat').transpose()
    bin_centres = dat[2]/1e3  # to Jy
    counts = dat[5]
    count_errors_l = dat[6]
    count_errors_u = dat[7]
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors_u = count_errors_u*(f/frq)**(1.5*alpha)
    count_errors_l = count_errors_l*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, [count_errors_l, count_errors_u]


def get_Hopkins_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/hopkins_2002.dat').transpose()
    bin_centres = dat[2]/1e3  # to Jy
    counts = dat[5]
    count_errors = dat[6]
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors

def get_Mauch_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/Mauch+2013.dat').transpose()
    bin_centres = dat[0]/1e3  # to Jy
    counts = dat[1]
    count_errors_lo = dat[2]-dat[1]
    count_errors_hi = dat[1]-dat[3]
    # flux Euc_count countlo counthi
    #counts = counts*(1e3)**2.5


    frq = 325.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors_lo = count_errors_lo*(f/frq)**(1.5*alpha)
    count_errors_hi = count_errors_hi*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, [count_errors_hi, count_errors_lo]


def get_Simpson_counts(f, alpha=-0.8):
    dat = np.loadtxt('/local/wwilliams/projects/src_counts/simpson_2006.dat').transpose()
    bin_centres = dat[2]/1e3  # to Jy
    counts = dat[5]
    count_errors = dat[6]
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors


def get_Bondi_counts(f, alpha=-0.8):
#Slow Shigh Smid N dNdS count Ecount C Nc
#0.0600 0.0735 0.066 380 9.24e10 3.32 0.17 0.99 1669 41
#0.0735 0.0900 0.081 324 6.43e10 3.84 0.21 1.16 1515 39
#0.0900 0.1103 0.100 203 3.29e10 3.26 0.23 1.36 1336 37
#0.1103 0.1351 0.122 142 1.88e10 3.09 0.26 1.60 1246 35
#0.1351 0.1655 0.150 135 1.46e10 3.99 0.34 1.29 822 29
#0.1655 0.2028 0.183 109 9.61e9 4.36 0.42 1.15 577 24
#0.2028 0.2484 0.224 74 5.33e9 4.02 0.47 1.00 393 20
#0.2484 0.3043 0.275 72 4.23e9 5.30 0.63 1.00 319 18
#0.3043 0.3728 0.337 45 2.16e9 4.49 0.67 1.00 247 16
#0.3728 0.4566 0.413 33 1.29e9 4.47 0.78 1.00 202 14
#0.4566 0.5593 0.505 28 8.95e8 5.14 0.97 1.00 169 13
#0.5593 0.6851 0.619 24 6.26e8 5.97 1.28 1.00 141 12
#0.6851 0.8393 0.758 12 2.56e8 4.05 1.17 1.00 117 11
#0.8393 1.0282 0.929 12 2.09e8 5.49 1.58 1.00 105 10

    bin_centres = np.array([0.066, 0.081, 0.1, 0.122, 0.15, 0.183, 0.224, 0.275, 0.337, 0.413, 0.505, 0.619, 0.758, 0.929])
    # mJy
    bin_centres = bin_centres/1e3  # to Jy
    counts = np.array([3.32, 3.84, 3.26, 3.09, 3.99, 4.36, 4.02, 5.3, 4.49, 4.47, 5.14, 5.97, 4.05, 5.49])
    count_errors = np.array([0.17, 0.21, 0.23, 0.26, 0.34, 0.42, 0.47, 0.63, 0.67, 0.78, 0.97, 1.28, 1.17, 1.58])
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors

def get_FIRST_counts(f, alpha=-0.8):
    # http://iopscience.iop.org/0004-637X/475/2/479/fulltext/
    bin_centres = np.array([0.001, 0.00112, 0.00126, 0.00141, 0.00158, 0.00178, 0.002, 0.00224, 0.00251, 0.00282, 0.00316, 0.00355, 0.00398, 0.00447, 0.00501, 0.00562, 0.00631, 0.00708, 0.00794, 0.00891, 0.01, 0.01122, 0.01259, 0.01413, 0.01585, 0.01778, 0.01995, 0.02239, 0.02512, 0.02818, 0.03162, 0.03548, 0.03981, 0.04467, 0.05012, 0.05623, 0.0631, 0.07079, 0.07943, 0.08913, 0.1, 0.1122, 0.12589, 0.14125, 0.15849, 0.17783, 0.19953, 0.22387, 0.25119, 0.28184, 0.31623, 0.35481, 0.39811, 0.44668, 0.50119, 0.56234, 0.63096, 0.70795, 0.79433, 0.89125, 1])
    # Jy
    counts = np.array([3.55, 3.23, 5.19, 7.23, 7.64, 9.12, 10.04, 12.27, 12.68, 14.16, 15.41, 17.1, 18.99, 20.54, 22.89, 25.26, 27.28, 31.38, 34.54, 37.24, 42.04, 46.01, 50.67, 54.59, 59.78, 66.18, 77.71, 85.57, 80.12, 100.88, 114.91, 102.91, 134.73, 119.97, 138.13, 160.17, 179.24, 172.49, 208.6, 230.73, 223.84, 231.56, 258.5, 246.72, 328.63, 260.16, 308.06, 321.1, 371.5, 375.95, 305.91, 296.63, 340.95, 351.93, 364.21, 328.27, 632.16, 482.54, 273.39, 490.07, 288.7])
    # Jy**1.5 / sr
    count_errors = np.array([0.04, 0.05, 0.06, 0.08, 0.09, 0.11, 0.12, 0.15, 0.17, 0.19, 0.22, 0.25, 0.29, 0.33, 0.37, 0.43, 0.49, 0.57, 0.65, 0.74, 0.85, 0.97, 1.11, 1.26, 1.43, 1.65, 1.94, 2.22, 2.35, 2.87, 3.34, 3.45, 4.3, 4.42, 5.17, 6.07, 7, 7.49, 8.98, 10.29, 11.05, 12.26, 14.12, 15.03, 18.92, 18.35, 21.77, 24.23, 28.41, 31.16, 30.64, 32.89, 38.44, 42.58, 47.22, 48.88, 73.94, 70.43, 57.79, 84.35, 74.38])
    # Jy**1.5 / sr

    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres, counts, count_errors


def get_gmrt_counts(f, alpha=-0.8, scale=False):
    [bin_centres, counts, count_errors_l, count_errors_u] = np.load('/local/wwilliams/phd/gmrt_bootes/bootes_analysis/fits/GMRT_BOOTES_153MHz_source_counts.npy')
    
    
    frq = 153.
    
    if scale:
        # scale to required freq
        bin_centres = bin_centres*(f/frq)**alpha
        counts = counts*(f/frq)**(1.5*alpha)
        count_errors_l = count_errors_l*(f/frq)**(1.5*alpha)
        count_errors_u = count_errors_u*(f/frq)**(1.5*alpha)
    
    
    return frq, bin_centres, counts, [count_errors_l, count_errors_u]
    

def  get_skads_counts_file(f, alpha=-0.8):
    data = pl.array(pl.loadtxt('/local/wwilliams/phd/gmrt_bootes/bootes_analysis/fits/skads_counts_model.dat')).transpose()
    bin_centres = pl.array(data[0])
    counts = pl.array(data[1])
    count_errors = pl.array(data[2])
    frq = 150.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
  
    return frq, bin_centres ,counts, count_errors


def smooth(x,window_len=11,window='hanning'):
  import numpy
  if x.ndim != 1:
      raise ValueError, "smooth only accepts 1 dimension arrays."

  if x.size < window_len:
      raise ValueError, "Input vector needs to be bigger than window size."

  if window_len<3:
      return x

  if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
      raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

  s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
  #print(len(s))
  if window == 'flat': #moving average
      w=numpy.ones(window_len,'d')
  else:
      w=eval('numpy.'+window+'(window_len)')

  y=numpy.convolve(w/w.sum(),s,mode='same')
  y=y[window_len/2+2:-window_len/2-1]
  return y

def get_areas(incat,cat_areafile,rmsfits,corFlux=1):
    print 'getting areas'
    #rmsdat = 0.79*pf.getdata(rmsfits)
    rmsdat = pf.getdata(rmsfits)
    rmsdat = rmsdat[np.isfinite(rmsdat)]
    rmsdat = rmsdat[rmsdat>0]
    rmsdat = rmsdat.flatten()
    rmsdat = rmsdat/corFlux
    rmshead = pf.getheader(rmsfits)
    ps = rmshead.get('CDELT2')

    areas = -1.*np.ones(len(incat))
    pixels = -1.*np.ones(len(incat))
    for i in range(len(incat)):
        fp = incat['Peak_flux'][i]
        rmslim = fp/5.
        pixel = pl.sum(rmsdat < rmslim)
        area = pixel*ps*ps
        areas[i] = area
        pixels[i] = pixel
    areas[pixels==0] = np.nan
    pixels[pixels==0] = np.nan

    incat = append_field(incat, 'Pixels', pixels)
    incat = append_field(incat, 'Area', areas)
    pf.writeto(cat_areafile,incat)
    return

def get_areas2(incat,cat_areafile,rmsfits,corFlux=1):
    print 'getting areas'
    #rmsdat = 0.79*pf.getdata(rmsfits)
    rmsdat = pf.getdata(rmsfits)
    rmsdat = np.ma.masked_where(np.isnan(rmsdat), rmsdat).compressed()
    rmsdat = rmsdat/corFlux
    rmshead = pf.getheader(rmsfits)
    ps = rmshead.get('CDELT2')

    areas = -1.*np.ones(len(incat))
    pixels = -1.*np.ones(len(incat))
    for i in range(len(incat)):
        fp = incat['Peak_flux'][i]
        rmslim = fp/5.
        pixel = pl.sum(rmsdat < rmslim)
        area = pixel*ps*ps
        areas[i] = area
        pixels[i] = pixel

    incat = append_field(incat, 'Pixels', rms)
    incat = append_field(incat, 'Pixels', pixels)
    incat = append_field(incat, 'Area', areas)
    pf.writeto(cat_areafile,incat)
    return


def plot_src_counts_all():
    
    fig, ax1 = pp.paper_double_ax()
    pl.minorticks_on()
    
    gmrt_frq, gmrt_bins, gmrt_counts, gmrt_count_errors = get_gmrt_counts(nu, alpha=alpha, scale=True)
    ax1.errorbar(gmrt_bins*scaleF, gmrt_counts, gmrt_count_errors, fmt='ko', color='gray', ms = 7, label='Williams et al. 2013')
    counts=counts_cor
    
    
    othercol = pl.cm.Set1(np.linspace(0, 1, 10))
    othersymb = ['x','*','o','^','v','>','<','d','p','*','h','H','1','2','3','4','D']
    othersymbi = 0
    FIRST_frq, FIRST_bins, FIRST_counts, FIRST_count_errors = get_FIRST_counts(nu, alpha=alpha)
    ax1.errorbar(FIRST_bins*scaleF, FIRST_counts, FIRST_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='White et al. 1997')
    othersymbi += 1
    
    Ciliegi_frq, Ciliegi_bins, Ciliegi_counts, Ciliegi_count_errors = get_Ciliegi_counts(nu, alpha=alpha)
    ax1.errorbar(Ciliegi_bins*scaleF, Ciliegi_counts, Ciliegi_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Ciliegi et al. 1999')
    othersymbi += 1
    
    
    
    Prandoni_frq, Prandoni_bins, Prandoni_counts, Prandoni_count_errors = get_Prandoni_counts(nu, alpha=alpha)
    ax1.errorbar(Prandoni_bins*scaleF, Prandoni_counts, Prandoni_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Prandoni et al. 2001')
    othersymbi += 1
    
    Richards_frq, Richards_bins, Richards_counts, Richards_count_errors = get_Richards_counts(nu, alpha=alpha)
    ax1.errorbar(Richards_bins*scaleF, Richards_counts, Richards_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Richards 2000')
    othersymbi += 1
    
    
    Hopkins_frq, Hopkins_bins, Hopkins_counts, Hopkins_count_errors = get_Hopkins_counts(nu, alpha=alpha)
    ax1.errorbar(Hopkins_bins*scaleF, Hopkins_counts, Hopkins_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Hopkins et al. 2002')
    othersymbi += 1
    
    
    Bondi_frq, Bondi_bins, Bondi_counts, Bondi_count_errors = get_Bondi_counts(nu, alpha=alpha)
    ax1.errorbar(Bondi_bins*scaleF, Bondi_counts, Bondi_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Bondi et al. 2003')
    othersymbi += 1
    
    Seymour_frq, Seymour_bins, Seymour_counts, Seymour_count_errors = get_Seymour_counts(nu, alpha=alpha)
    ax1.errorbar(Seymour_bins*scaleF, Seymour_counts, Seymour_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Seymour et al. 2004')
    othersymbi += 1
    
    Simpson_frq, Simpson_bins, Simpson_counts, Simpson_count_errors = get_Simpson_counts(nu, alpha=alpha)
    ax1.errorbar(Simpson_bins*scaleF, Simpson_counts, Simpson_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Simpson et al. 2006')
    othersymbi += 1
    
    Mauch_frq, Mauch_bins, Mauch_counts, Mauch_count_errors = get_Mauch_counts(nu, alpha=alpha)
    ax1.errorbar(Mauch_bins*scaleF, Mauch_counts, Mauch_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='Mauch et al. 2013')
    othersymbi += 1
    
    
    w_frq, w_bins, w_counts, w_count_errors = get_deZotti_counts(nu, alpha=alpha, value='wi97')
    ax1.errorbar(w_bins*scaleF, w_counts, w_count_errors, color='purple', fmt='o', ms = 5,label='zotti w')
    
    w_frq, w_bins, w_counts, w_count_errors = get_deZotti_counts(nu, alpha=alpha, value='all')
    ax1.errorbar(w_bins*scaleF, w_counts, w_count_errors, color='gray', fmt='.', ms = 5,label='zotti all')
    
    
    y1,y2 = ax1.get_ylim()
    
    
    #pl.title(' Source Counts ')
    pl.ylabel(r'$S^{5/2} dN / dS$ [Jy$^{3/2}$ sr$^{-1}$]')
    pl.xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='medium'
    pl.legend(loc=4,numpoints=1,frameon=False)
    #pl.minorticks_on()
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))
    
    ax1.set_xlim(0.15,9000)
    ax1.set_ylim(2, 9000)
    
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    
    pp.fig_save_many(fig, '%s_source_counts_comp' %(name))
    
    return

#def S14(S150):
    ## S150 in Jy, S14 in mJy
    #return 1000*S150*(1400./150)**-0.8

#def phimed(S14):
    ## S14 in mJy
    #return 2*(S14)**0.3

#def h(phi,S150):                  
    #S = S14(S150)
    ##return 2**(-1*(phi/phimed(S))**0.62)
    #return np.exp(np.log(2)*(phi/phimed(S))**0.62)


#def h(phi,S,th=1):    
    ##return 2**(-1*(phi/phimed(S))**0.62)
    #return np.exp(-np.log(2)*(phi/(th*phimed(S)))**0.62)


#def h2(phi,S,th=1):    
    ##return 2**(-1*(phi/phimed(S))**0.62)
    #return np.exp(-np.log(2)*(phi/(phimed(S)))**(0.62))

##def h2(theta):
    ##f = (1./1.6**theta)*(theta<4) +  (theta**-1.3 -0.01)*(theta>=4)
    ##f[f<0] = 0
    ##return f
  
def S14(S150, alpha=-0.8):
    # S150 in Jy, S14 in mJy
    return 1000*S150*(1400./150)**alpha

def phimed(S14):
    # S14 in mJy
    return 2*(S14)**0.3

def phimed2(S14):
    # S14 in mJy
    return 2.*(S14 < 1.) +  (2*(S14)**0.3)*(S14>=1.)

def h(phi,S150,th=1):                  
    S = S14(S150)
    #return 2**(-1*(phi/phimed(S))**0.62)
    return np.exp(-1*np.log(2)*(phi/phimed(S))**0.62)

def h2(phi,S150):                  
    S = S14(S150)
    #return 2**(-1*(phi/phimed(S))**0.62)
    return np.exp(-1*np.log(2)*(phi/phimed2(S))**0.62)

def dh2(phi,S150):                  
    S = S14(S150)
    #return 2**(-1*(phi/phimed(S))**0.62)
    return np.exp(-1*np.log(2)*(phi/phimed2(S))**0.62) * (0.62*np.log(2)/phimed2(S)) * (phi/phimed2(S))**-1.62
  
  
  
def get_res_cor(S):
    srange, RC1low, RC1high, RC2low, RC2high = np.load('rescor.npy')
    si = sum(srange<S)
    if si == len(srange):
        si -= 1
    return np.array([RC1low[si], RC1high[si], RC2low[si], RC2high[si]])

  
def get_comp_cor(S):
    area_i,irangep,C,C_sig,C_area = np.load('LOFAR_BOOTES_Detection_fraction_int_i.npy')
    irangep /= 1000.
    si = sum(irangep<S)
    if si == len(irangep):
        si -= 1
    return np.array(1./C_area[si])

  
  
def get_fdr_cor(S):
    irangep,C,C_sig = np.load('LOFAR_BOOTES_Rely_fraction_1sig_i.npy')
    irangep /= 1000.
    si = sum(irangep<S)
    if si == len(irangep):
        si -= 1
    Ci = C[si]
    if np.isnan(Ci):
        Ci = 0.
    return np.array(1.- Ci)

  
#def plot_source_counts(sourcecat, rmsfits, name, nu, corFlux=1., unit='mJy', clobber_areas=False):
    
    
if 1:
    
    parser = argparse.ArgumentParser() 
    parser.add_argument('fits_cat', help="Name of fits srouce catalogue")
    parser.add_argument('rmsfits', help="Name of rms fits image for area")
    parser.add_argument('name', help="Identifying Name")
    parser.add_argument('nu', help="Frequency")
    parser.add_argument('--unit', default='mJy', help="(default mJy)")
    parser.add_argument('-c','--clobber', action='store_true', default=False, help="(default false)")
    parser.add_argument('-s','--scalepeak', action='store_true', default=False, help="(default false)")
    parser.add_argument('-m','--multiple_correction', action='store_true', default=False, help="(default false)")
    
    args = parser.parse_args()
    
    fits_cat = args.fits_cat
    rmsfits = args.rmsfits
    name = args.name
    nu = float(args.nu)
    corFlux=1. 
    unit=args.unit
    #clobber_areas=args.clobber
    sourcecat = fits_cat
    
    #clobber_areas=False
    
    
    alpha = -0.8
    if unit=='mJy':
        scaleF = 1000.
    elif unit=='Jy':
        scaleF = 1.
    else:
        print 'ERROR: unit not handled: '+unit
    
    cat = load_fits(sourcecat )
    #rmsfits =  'A2256.MOSAIC.MASKED.RMS.FITS'

    
    #cat_areafile =name+'.areas.fits'
    #if not os.path.isfile(cat_areafile):
        #get_areas(cat,cat_areafile,rmsfits)
    #elif clobber_areas:
        #os.system('rm '+cat_areafile)
        #get_areas(cat,cat_areafile,rmsfits)
    #cat_areas = load_fits(cat_areafile )

    if "Flag_artefact" in cat.dtype.names:
        nbad = np.sum(np.asarray(cat["Flag_artefact"])!=0)
        cat = cat[cat["Flag_artefact"]==0]
        print "deselecting {n} sources flagged as artefacts".format(n=nbad)



    
    if args.multiple_correction:
        
        print "doing multiple source correction - v0"
        
        matchcoord=ac.SkyCoord(cat['RA'],cat['DEC'],unit='deg')
        catcoord=ac.SkyCoord(cat['RA'],cat['DEC'],unit='deg')
        idx,sep,d = ac.match_coordinates_sky(matchcoord, catcoord, nthneighbor=2)
        fluxrat = cat['Total_flux']/cat['Total_flux'][idx]
        fluxsum = cat['Total_flux'] + cat['Total_flux'][idx]
        
        stheta = 1e-3*(1400./150)**0.8*10*(sep.arcsec/100)**2.  # in Jy # magliocchetti+ 1998
        tt2 = (stheta < fluxsum)

        tt = (fluxrat>= 0.25) & (fluxrat <= 4)

        idx0 = np.arange(len(cat))

        idxm = idx[tt&tt2]
        idx0m = idx0[tt&tt2]
        print "{n} match criteria to be joined (double counted)".format(n=len(idxm))
        
        # replace flux with sum for multiples
        cat['Total_flux'][idxm] = fluxsum[idxm]
        
        # now make sure we remove the multiples - note this works for pairs...
        selm = []
        for i in range(len(cat)):
            if i in idxm:
                mi = np.where(idxm==i)[0]
                i2 = idx0m[mi[0]]
                #print i, mi, i2
                if (i not in selm) and (i2 not in selm):
                    selm.append(i)
            else:
                selm.append(i)
                
        print "{n} sources".format(n=len(cat))
        print "{n} selected".format(n=len(selm))
        cat = cat[selm]


    RA0 = 218.
    DEC0 = 34.5
    rad = angular_separation(cat['RA'], cat['DEC'],RA0,DEC0)
    bloss = bandwidth_smearing2(97.5e3,150.e6,7.5,rad*3600.)
    tloss = time_smearing2(8, 5.5,rad*3600)
    loss = bloss*tloss
    #loss2 = 0.833
    loss2 = 1.
    peakcor = cat['Peak_flux'] / (loss*loss2)

    if not args.scalepeak:
        PSNRLIM = 5
        PSNRLIM_DET = 5
        FSNRLIM = 5.
        sel = cat['Peak_flux']/cat['Isl_rms'] > PSNRLIM
        
        nsel = np.sum(sel)
        nall = len(sel)
        
        print "{n1} of {n2} selected on apparent peak flux < {plim}".format(n1=nsel, n2=nall, plim=PSNRLIM)
    
    else:
        print "scaling peak"
        PSNRLIM_DET = 5.0*1.2
        PSNRLIM = 5.
        
        PSNRLIM_DET = 5.0
        PSNRLIM = 5.*1.2
        
        PSNRLIM_DET = 5.0
        PSNRLIM = 5.*1.2
        
        FSNRLIM = 5.
        
        sel = peakcor/cat['Isl_rms'] > PSNRLIM
        selreal = cat['Peak_flux']/cat['Isl_rms'] > PSNRLIM
    
        nsel = np.sum(sel)
        nselreal = np.sum(selreal)
        nall = len(sel)
        print "{n1} of {n2} selected on apparent peak flux < {plim}".format(n1=nselreal, n2=nall, plim=PSNRLIM)
        print "{n1} of {n2} selected on 'corrected' peak flux < {plim}".format(n1=nsel, n2=nall, plim=PSNRLIM)
        
    peakcor = peakcor[sel]
    
    bloss = bloss[sel]
    tloss = tloss[sel]
    
    #cat['DC_Maj'] = cat['DC_Maj'] * bloss  # correct sizes
    #cat['DC_Min'] = cat['DC_Min'] * tloss  
    
    cat = cat[sel]
    

    
    print '# source counts #'
    print 'psnrlim_det:',PSNRLIM_DET
    print 'psnrlim:',PSNRLIM
    print 'fsnrlim:',FSNRLIM
    
    #Area = 30.   # square degrees
    
    # HARDCODED #
    app=(1.5/3600)**2.  # area per pixel in sq deg
    degtosterad = (pl.pi/180.)**2.
    #Area = Area * degtosterad
    
    #pixel scale : pix = ps degrees

    ra_cen, dec_cen = 218.02396, 34.279861
    # Make plot of extrapolated/measured flux vs. distance from phase center
    
    #flux =  np.asarray(cat['Total_flux'].copy())
    flux =  np.asarray(cat['Isl_Total_flux'].copy())
    peak =  np.asarray(cat['Peak_flux'].copy())
    noise =  np.asarray(cat['Isl_rms'].copy())
    eflux =  np.asarray(cat['E_Total_flux'].copy())
    #areap = np.asarray(cat_areas['Pixels'].copy())
    
    
    
    # correct poisson errors
    PE = {"S": 1, "beta": 0., "gamma": 1.}
    #rmslim = bins[-1]/5.
    
    
    if 'fits' in rmsfits:
        import pyfits as pf
        
        print "get rmsmap"
        rmsdat = pf.getdata(rmsfits)
        #rmsdat = np.ma.masked_where(np.isnan(rmsdat), rmsdat).compressed()
        rmsdat = rmsdat[np.isfinite(rmsdat)]
        #rmsdat = rmsdat /corFlux
        rmshead = pf.getheader(rmsfits)
        ps = rmshead.get('CDELT2')
        
        rmspixels = rmsdat[rmsdat>0]

        pixels_tot =  len(rmspixels) # pl.sum(rmsdat < rmslim)
        Area_tot = pixels_tot*ps*ps
    
        print "get areas"
        #n,b=np.histogram(rmspixels,  bins=1000, range=(0,0.001))
        n,logb=np.histogram(np.log10(rmspixels),  bins=100, range=(np.log10(90e-6), np.log10(350e-6)))
        logbc = (10**logb[1:] + 10**logb[:-1])/2.
        nc = np.cumsum(n)
        nc = 1.*nc/nc[-1]  # normalise
        rms90 = logbc[np.sum(nc<0.9)]
        rms50 = logbc[np.sum(nc<0.5)]
        #rms20 = b[np.sum(nc<0.2)]
        rms10 = logbc[np.sum(nc<0.1)]
        print "log version"
        print rms10, rms50, rms90 
        
        #n,b=np.histogram(rmspixels,  bins=500, range=(100.e-6,0.001))
        rmspixels1=rmspixels[(rmspixels<0.001)&(rmspixels>100.e-6)]
        n,b=np.histogram(rmspixels1,  bins=500)
        bc = (b[1:] + b[:-1])/2.
        nc = np.cumsum(n)
        nc = 1.*nc/nc[-1]  # normalise
        
        #rms90 = bc[np.sum(nc<0.9)]
        ##rms75 = bc[np.sum(nc<0.75)]
        #rms50 = bc[np.sum(nc<0.5)]
        rms25 = b[np.sum(nc<0.25)]
        #rms10 = bc[np.sum(nc<0.1)]
        print "linear version"
        #print rms10, rms50, rms90 
        
        noiserange = np.logspace(np.log10(90e-6), np.log10(350e-6), 25)
        areas = np.ones_like(noiserange)
        for i,rms in enumerate(noiserange):
            areas[i] = np.sum(rmspixels<rms)*app
        
        
        import pickle
        rmsdict = {"b": b, "nc": nc}
        with open('rms.pickle','w') as picklefile:
            pickle.dump(rmsdict,picklefile)
        

    elif 'pickle' in rmsfits:
        import pickle
        rmsdict = pickle.load(rmsfits)
        b = rmsdict['b']
        nc = rmsdict['nc']



    if not args.scalepeak:
        #areaf_av=np.array([ nc[np.min((np.sum(b<ip/(5.*1.4)),len(nc)-1))] for ip in peak])
        areaf = np.array([ nc[np.min((np.sum(b<ip/(PSNRLIM)),len(nc)-1))] for ip in peak])
        areafscale = np.array([ nc[np.min((np.sum(b<ip/(PSNRLIM_DET)),len(nc)-1))] for ip in peak])
        #areaf = areap/pixels_tot
    else:
        areaf = np.array([ nc[np.min((np.sum(b<ip/(PSNRLIM)),len(nc)-1))] for ip in peakcor])
        areafscale = np.array([ nc[np.min((np.sum(b<ip/(PSNRLIM_DET)),len(nc)-1))] for ip in peak])
        

    print "get bins"

    #b1 = PSNRLIM_DET*rms90
    #b1 = PSNRLIM_DET*rms75
    b1 = PSNRLIM_DET*rms25
    b2 = 0.25*flux.max()
    #bmid = 10**((pl.log10(b1) + pl.log10(b2))/2.)
    #nb1 = 30
    #nb2 = 10
    #nb=50
    ##nb3 = 10
    #bins1 = pl.logspace(pl.log10(b1),pl.log10(bmid),nb1)
    #bins2 = pl.logspace(pl.log10(bmid),pl.log10(b2),nb2)
    #bins = np.hstack((bins1,bins2))
    #bins = np.unique(bins)
    #bins = pl.logspace(pl.log10(b1),pl.log10(b2),nb)
    
    A = 51.18
    beta = -1.436
    #beta = -1.6
    Nperbin = 200.
    Nperbin = 150.
    bins=[b1]
    while bins[-1] < b2: 
        bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
    
    
    
    nbins = len(bins) - 1 
    bin_centres = pl.zeros( nbins )
    bin_widths = pl.zeros( nbins )
    bin_errl = pl.zeros( nbins )
    bin_erru = pl.zeros( nbins )
    for i in range(nbins):
        bin_centres[ i ] =  np.sqrt(bins[i]*bins[i+1])  #(bins[i] + bins[i+1])/2.
        bin_widths[ i ] = (bins[i+1] - bins[i])
        bin_erru[ i ] = (bins[i+1] - bin_centres[ i ])
        bin_errl[ i ] = (bin_centres[ i ] - bins[i])
    bincount = pl.zeros( nbins )
    bincountpeak = pl.zeros( nbins )
    #dnbins = pl.zeros( nbins )
    counts_areascale = pl.zeros( nbins )
    counts = pl.zeros( nbins )
    counts_low = pl.zeros( nbins )
    counts_high = pl.zeros( nbins )
    counts_cor = pl.zeros( nbins )
    counts_weight = pl.zeros( nbins )
    counts_weight1 = pl.zeros( nbins )
    counts_weight2 = pl.zeros( nbins )
    counts_area = pl.zeros( nbins )
    count_errors = pl.zeros( nbins )
    count_errors_l = pl.zeros( nbins )
    count_errors_u = pl.zeros( nbins )
    count_cor_errors = pl.zeros( nbins )
    corF = pl.ones(nbins)
    corFe = pl.ones(nbins)
    counts_corcomp = pl.zeros( nbins )
    counts_corfdr = pl.zeros( nbins )
    counts_corres = pl.zeros( nbins )
    counts_corres_err = pl.zeros( nbins )
    counts_correslow = pl.zeros( nbins )
    counts_correshigh = pl.zeros( nbins )

    print bins
    print nbins
    print len(bins)
    print 'Bin [mJy]   & <S> & <W> & Raw Count + - & Cor Count & Area L C R & Norm Count | L R & Norm Count Errors'

    
    
    ################### RAW counts ###
    
    for i in range( nbins ): 
        mask = (flux>bins[i]) * (flux<=bins[i+1]) *(flux/noise >= FSNRLIM ) #*(peak/noise >= PSNRLIM )
        cnt = pl.sum( mask )
        bincount[ i ] = cnt
    for i in range( nbins ): 
        mask = (peak>bins[i]) * (peak<=bins[i+1]) #*(peak/noise >= PSNRLIM )
        cnt = pl.sum( mask )
        bincountpeak[ i ] = cnt
        
        
    print "plot raw counts"
    
    area_per_bin = np.zeros(nbins)
    area_per_bin2 = np.zeros(nbins)
    for i,bini in enumerate(bin_centres):
        area_per_bin[i] = np.sum(rmspixels<(bini/PSNRLIM_DET))*app
    area_per_bin = area_per_bin/area_per_bin[-1]
        
        
    for i,bini in enumerate(bin_centres):
        ii = np.min((np.sum(bc<(bini/PSNRLIM_DET)),len(nc)-1))
        area_per_bin2[i] =  nc[ii]
        
    fig, ax1 = pp.paper_single_ax()
    ax1.semilogx(bin_centres*scaleF, bincount,  'ko', ms = 10, zorder=10)
    ax1.errorbar(bin_centres*scaleF, bincount, np.sqrt(bincount), [bin_erru*scaleF,bin_errl*scaleF],  'ko', ms = 10, label='This work', zorder=10)
    ax2 = ax1.twinx()
    ax2.plot(bin_centres*scaleF, area_per_bin,'b')
    ax2.plot(bin_centres*scaleF, area_per_bin,'r')
    ax2.set_ylim(0,1.01)
    ax2.set_ylabel(r'Area')
    ax1.set_ylabel(r'N')
    ax1.set_xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='medium'
    pl.legend(loc=4,numpoints=1,frameon=False)
    pl.minorticks_on()
    pl.subplots_adjust(right=0.8)
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    pp.fig_save_many(fig, '%s_source_counts_raw' %(name))
    fig, ax1 = pp.paper_single_ax()
    
    fig, ax1 = pp.paper_single_ax()
    ax1.loglog(bin_centres*scaleF, (1/area_per_bin)*bincount/bin_widths,  'ko', ms = 10)
    ax1.errorbar(bin_centres*scaleF, (1/area_per_bin)*bincount/bin_widths, (1/area_per_bin)*np.sqrt(bincount)/bin_widths, [bin_erru*scaleF,bin_errl*scaleF],  'ko', ms = 10, label='This work')
    ax2 = ax1.twinx()
    ax2.plot(bin_centres*scaleF, area_per_bin,'b')
    ax2.plot(bin_centres*scaleF, area_per_bin,'r')
    ax2.set_ylim(0,1.01)
    ax2.set_ylabel(r'Area')
    ax1.set_ylabel(r'$dN/dS$')
    ax1.set_xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='medium'
    pl.legend(loc=4,numpoints=1,frameon=False)
    pl.minorticks_on()
    pl.subplots_adjust(right=0.8)
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    pp.fig_save_many(fig, '%s_source_counts_dnds_raw' %(name))
    
    
    #gmrt_bins
    
    b1 = 1.5e-2
    b2 = 3.
    nb = 20
    gmrt_bins = pl.zeros(nb)
    gmrt_bins[0] = b1
    for i in pl.arange(1,nb):
        bf = 4.**0.125
        if gmrt_bins[i-1] > 0.086:
            bf = 4.**0.25
        if gmrt_bins[i-1] > 0.6:
            bf = 4.**0.5
        if gmrt_bins[i-1] > 1.0:
            bf = 4.
        gmrt_bins[i] = gmrt_bins[i-1]*bf
    gmrt_bin_centres = pl.zeros(nb-1)
    gmrt_bin_widths = pl.zeros(nb-1)
    gmrt_bincount = pl.zeros(nb-1)
    for i in range(nb - 1):
        gmrt_bin_centres[i] = (gmrt_bins[i] + gmrt_bins[i+1])/2.
        gmrt_bin_widths[i] = gmrt_bins[i+1] - gmrt_bins[i]
    
    for i in range( nb-1 ): 
        mask = (flux>gmrt_bins[i]) * (flux<=gmrt_bins[i+1]) *(flux/noise >= FSNRLIM )*(peak/noise >= PSNRLIM )
        cnt = pl.sum( mask )
        gmrt_bincount[ i ] = cnt
    gmrt_area_per_bin = np.zeros(nb-1)
    for i,bini in enumerate(gmrt_bin_centres):
        gmrt_area_per_bin[i] = np.sum(rmspixels<(bini/PSNRLIM_DET))*app
    gmrt_area_per_bin = gmrt_area_per_bin * degtosterad
    
    fig, ax1 = pp.paper_single_ax()
    ax1.loglog(gmrt_bin_centres*scaleF, gmrt_bin_centres**2.5*(1/gmrt_area_per_bin)*gmrt_bincount/gmrt_bin_widths,  'ko', ms = 10)
    ax1.errorbar(gmrt_bin_centres*scaleF, gmrt_bin_centres**2.5*(1/gmrt_area_per_bin)*gmrt_bincount/gmrt_bin_widths, gmrt_bin_centres**2.5*(1/gmrt_area_per_bin)*np.sqrt(gmrt_bincount)/gmrt_bin_widths,  fmt='ko', ms = 10, label='This work')
    ax1.set_ylabel(r'$dN/dS$')
    ax1.set_xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='medium'
    pl.legend(loc=4,numpoints=1,frameon=False)
    pl.minorticks_on()
    pl.subplots_adjust(right=0.8)
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    pp.fig_save_many(fig, '%s_source_counts_gmrt_bins' %(name))
    save_gmrt_bins = gmrt_bin_centres
    save_gmrt_counts = gmrt_bin_centres**2.5*(1/gmrt_area_per_bin)*gmrt_bincount/gmrt_bin_widths
    save_gmrt_count_errs = gmrt_bin_centres**2.5*(1/gmrt_area_per_bin)*np.sqrt(gmrt_bincount)/gmrt_bin_widths
    
    
    fig, ax1 = pp.paper_single_ax()
    ax1.semilogx(bin_centres*scaleF, bincountpeak,  'ko', ms = 10)
    ax1.errorbar(bin_centres*scaleF, bincountpeak, np.sqrt(bincount), [bin_erru*scaleF,bin_errl*scaleF],  'ko', ms = 10, label='This work')
    ax2 = ax1.twinx()
    ax2.plot(bin_centres*scaleF, area_per_bin,'b')
    ax2.set_ylim(0,1.01)
    ax2.set_ylabel(r'Area')
    ax1.set_ylabel(r'N')
    ax1.set_xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='medium'
    pl.legend(loc=4,numpoints=1,frameon=False)
    pl.minorticks_on()
    pl.subplots_adjust(right=0.8)
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    pp.fig_save_many(fig, '%s_source_counts_peak_raw' %(name))
    
    print "done"
    
    #######################
    
    
    wdist_bins = 1e-3*np.array([5., 10.,20,30,40])
    for i in range( len(wdist_bins) -1 ): 
        mask = (flux>wdist_bins[i]) * (flux<=wdist_bins[i+1])
        noiselim1 = 0.11e-3
        noiselim2 = 0.15e-3
        BMAJ = 7.4/3600.
        BMIN = 5.6/3600.
        limA1 = 3600*np.sqrt(BMAJ*BMIN*wdist_bins[i]/(5*noiselim1))
        limA2 = 3600*np.sqrt(BMAJ*BMIN*wdist_bins[i]/(5*noiselim2))
        meanbeam = 3600.*np.sqrt(BMAJ*BMIN)
        DC_limA1 = np.sqrt(limA1**2. - meanbeam**2.)
        DC_limA2 = np.sqrt(limA2**2. - meanbeam**2.)
        DC_lim = np.min((DC_limA1,DC_limA2))
        ##limA = np.max((lim1,lim2,lim3), axis=0)
        ##ax.fill_between(srange*Fscale, limA1, limA2, color='b', alpha=0.2)
        #ax.fill_between(srange*Fscale, DC_limA1, DC_limA2, color='b', alpha=0.2)
        
        bmins = pl.ma.masked_where(1-mask,np.sqrt(cat['Min']*cat['Min'])*3600.).compressed()
        bmajs = pl.ma.masked_where(1-mask,np.sqrt(cat['Maj']*cat['Maj'])*3600.).compressed()
        bths = pl.ma.masked_where(1-mask,np.sqrt(cat['Maj']*cat['Min'])*3600.).compressed()
        bmajs = pl.ma.masked_where(1-mask,np.sqrt(cat['DC_Maj']*cat['DC_Maj'])*3600.).compressed()
        fig, axi = pp.paper_single_ax()
        axi.set_title('${b1:.2f} <S< {b2:.2f}$ mJy'.format(b1=wdist_bins[i]*scaleF, b2=wdist_bins[i+1]*scaleF))
        axi.vlines(3600.*BMAJ,0,1,colors='k', linestyles='dotted')
        axi.vlines(DC_lim,0,1,colors='k', linestyles='dotted')
        #axi.vlines(DC_limA2,0,1,colors='k', linestyles='dotted')
        #axi.hist(bmajs, normed=True, cumulative=-1, histtype='step', color='g', bins=np.linspace(0,100,101))
        #axi.hist(bmajs, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        #axi.hist(bmins, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        #axi.hist(bths, normed=True, cumulative=-1, histtype='step', color='k', bins=np.linspace(0,100,101))
        #bmajscor = pl.ma.masked_where(1-mask,np.sqrt(bloss*tloss*cat['DC_Maj']*cat['DC_Min'])*3600.).compressed()
        #axi.hist(bmajscor, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        #axi.hist(bmajscor, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        x1,x2 = axi.get_xlim()
        x1 = 0.
        x2 = 40.
        xr = np.linspace(x1,x2,101)
        hr1 = h(xr,wdist_bins[i])
        hr2 = h(xr,wdist_bins[i+1])
        #hr21 = h2(xr,wdist_bins[i])
        #hr22 = h2(xr,wdist_bins[i+1])
        axi.hist(bmajs, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        hr21 = h2(xr,wdist_bins[i])
        hr22 = h2(xr,wdist_bins[i+1])
        #axi.fill_between(xr,hr1,hr2,color='r',alpha=0.5)
        dxr = (xr[1:]-xr[:-1])
        N21 = np.sum(hr21[1:]/dxr)
        N22 = np.sum(hr22[1:]/dxr)
        N21 = 1
        N22 = 1

        axi.fill_between(xr,hr21/N21,hr22/N22,color='m',alpha=0.5)
        axi.set_xlabel(r'$\theta $ [arcsec]')
        axi.set_ylabel(r'$h(\theta>\theta_{\mathrm{max}}) $ [arcsec]')
        axi.set_xlim(0,40)
        axi.set_ylim(0,1)
        pp.fig_save_many(fig, '%s_maj_wbin%i' %(name,i))

    
        
    for i in range( nbins ): 
        mask = (flux>bins[i]) * (flux<=bins[i+1]) #*(peak/noise >= PSNRLIM )
        
        
        
        
        #*(flux/noise >= FSNRLIM )
        cnt = pl.sum( mask )
        cntw = pl.sum( mask )
        #fluxes = pl.ma.masked_where(1-mask,flux).compressed()
        #areaps = pl.ma.masked_where(1-mask,areap).compressed()
        #areaf_avs = pl.ma.masked_where(1-mask,areaf_av).compressed()
        areafscales = pl.ma.masked_where(1-mask,areafscale).compressed()
        areafs = pl.ma.masked_where(1-mask,areaf).compressed()
        #mean_f = pl.average(areaps)
        noises = pl.ma.masked_where(1-mask,noise).compressed()
        peaks = pl.ma.masked_where(1-mask,peak).compressed()
        fluxs = pl.ma.masked_where(1-mask,flux).compressed()
        cnt_err = pl.sqrt(cnt)
        # poisson errors ! low counts
        if cnt > 0:
            cnt_u_err = (cnt + PE['S']*pl.sqrt(cnt+0.75) + (PE['S']**2+3)/4.) - cnt
            cnt_l_err = cnt - ( cnt*(1. - 1./(9.*cnt) - PE['S']/(3.*pl.sqrt(cnt)) + PE['beta']*cnt**PE['gamma'] )**3. )
        else:
            cnt_u_err = 0
            cnt_l_err = 0
        # gaussian errors
        #cnt_u_err = pl.sqrt(cnt)
        #cnt_l_err = pl.sqrt(cnt)
        
        #npt, pt = np.histogram(peaks/fluxs, range=(0,1), bins=10)#,normed=True)
        #theta = np.sqrt(7.5*5.5/pt[1:])
        #fraclarge = h(theta, bins[i])
        #ncorpt = npt *(1+ fraclarge)
        
        #ptcrit = PSNRLIM_DET*np.array((rmspixels.min(),rmspixels.max())) / bin_centres[i]
        #ptcrit = PSNRLIM_DET*np.array((rms10,rms90)) / bin_centres[i]
        #thetacrit = np.sqrt(7.5*5.5/ptcrit)
        #fraclarge = h2(thetacrit)
        #COR_RES = np.sum(ncorpt)/np.sum(npt)
        #COR_RES = 1+fraclarge
        #COR_RES = np.max(COR_RES)
        COR_RES_arr = get_res_cor(bin_centres[i])
        
        
        # completeness correction
        counts_corcomp[i] = get_comp_cor(bin_centres[i])
        #print bin_centres[i], counts_corcomp[i]
        
        # fdr correction
        counts_corfdr[i] = get_fdr_cor(bin_centres[i])

        #COR_RES = 1
        #print COR_RES_arr
        COR_RES = COR_RES_arr[2]
        minCOR_RES = COR_RES_arr.min()
        maxCOR_RES = COR_RES_arr.max()
        
        counts_corres[i] = COR_RES
        counts_correslow[i] = minCOR_RES
        counts_correshigh[i] = maxCOR_RES
        counts_corres_err[i] = np.sqrt((0.1*COR_RES)**2.+ (maxCOR_RES-minCOR_RES)**2.)
        
        noiselim1 = 0.11e-3
        noiselim2 = 0.15e-3
        BMAJ = 7.4/3600.
        BMIN = 5.6/3600.
        limA1 = 3600*np.sqrt(BMAJ*BMIN*bins[i]/(5*noiselim1))
        limA2 = 3600*np.sqrt(BMAJ*BMIN*bins[i]/(5*noiselim2))
        meanbeam = 3600.*np.sqrt(BMAJ*BMIN)
        DC_limA1 = np.sqrt(limA1**2. - meanbeam**2.)
        DC_limA2 = np.sqrt(limA2**2. - meanbeam**2.)
        DC_lim = np.min((DC_limA1,DC_limA2))
        ##limA = np.max((lim1,lim2,lim3), axis=0)
        ##ax.fill_between(srange*Fscale, limA1, limA2, color='b', alpha=0.2)
        #ax.fill_between(srange*Fscale, DC_limA1, DC_limA2, color='b', alpha=0.2)
        
        bmajs = pl.ma.masked_where(1-mask,np.sqrt(cat['Maj']*cat['Min'])*3600.).compressed()
        bmajs = pl.ma.masked_where(1-mask,np.sqrt(cat['DC_Maj']*cat['DC_Min'])*3600.).compressed()
        fig, axi = pp.paper_single_ax()
        axi.set_title('${b1:.2f} <S< {b2:.2f}$ mJy'.format(b1=bins[i]*scaleF, b2=bins[i+1]*scaleF))
        axi.vlines(3600.*BMAJ,0,1,colors='k', linestyles='dotted')
        axi.vlines(DC_lim,0,1,colors='k', linestyles='dotted')
        #axi.vlines(DC_limA2,0,1,colors='k', linestyles='dotted')
        #axi.hist(bmajs, normed=True, cumulative=-1, histtype='step', color='g', bins=np.linspace(0,100,101))
        axi.hist(bmajs, normed=True, cumulative=-1, histtype='step', color='g', bins=np.linspace(0,100,101))
        
        bmajscor = pl.ma.masked_where(1-mask,np.sqrt(bloss*tloss*cat['DC_Maj']*cat['DC_Maj'])*3600.).compressed()
        #axi.hist(bmajscor, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        axi.hist(bmajscor, normed=True, cumulative=-1, histtype='step', color='b', bins=np.linspace(0,100,101))
        x1,x2 = axi.get_xlim()
        xr = np.linspace(x1,x2,101)
        hr1 = h(xr,bins[i])
        hr2 = h(xr,bins[i+1])
        hr21 = h2(xr,bins[i])
        hr22 = h2(xr,bins[i+1])
        axi.fill_between(xr,hr1,hr2,color='r',alpha=0.5)
        axi.fill_between(xr,hr21,hr22,color='m',alpha=0.5)
        axi.set_xlabel(r'$\theta $ [arcsec]')
        axi.set_xlim(0,40)
        axi.set_ylim(0,1)
        pp.fig_save_many(fig, '%s_maj_bin%i' %(name,i))


        
        # area of map in which source could be found (S > 5 sigma, sigma < S/5)
        rmslimc_scale = bin_centres[i]/PSNRLIM
        rmslimc = bin_centres[i]/PSNRLIM_DET
        rmsliml = bins[i]/PSNRLIM_DET
        rmslimr = bins[i+1]/PSNRLIM_DET
        pixelsc_scale = pl.sum(rmspixels < rmslimc_scale)
        pixelsc = pl.sum(rmspixels < rmslimc)
        pixelsl = pl.sum(rmspixels < rmsliml)
        pixelsr = pl.sum(rmspixels < rmslimr)
        #if pixelsc == 0:
            #continue
        Area_tot = pixels_tot*ps*ps
        Areac_areascale = pixelsc_scale*ps*ps
        Areac = pixelsc*ps*ps
        Areal = pixelsl*ps*ps
        Arear = pixelsr*ps*ps
        Area_stc_areasclae = Areac_areascale * degtosterad
        Area_stc = Areac * degtosterad
        Area_stl = Areal * degtosterad
        Area_str = Arear * degtosterad
        Area_sttot = Area_tot * degtosterad
        #if pixelsl < pixels:
        #if 1:
        #weights = pl.zeros(len(areaps))
        #weightsl = pl.zeros(len(areaps))
        #weightsr = pl.zeros(len(areaps))
        #weights = areaps/float(pixels_tot)
        #weightsl = areaps/float(pixelsl)
        #weightsr = areaps/float(pixelsr)
        #weightsc = areaps/float(pixelsc)
            #print f, rmslim, pixels, weights[fi]
        #cntc = pl.sum( Area_stc/(areaps/float(pixelsc)) )
        #cntr = pl.sum( Area_str/(areaps/float(pixelsr)) )
        #cntl = pl.sum( Area_stl/(areaps/float(pixelsl)) )
        #cntw = pl.sum( Area_stc/weights ) #/pl.sum()
        #cntl = pl.sum( 1./weightsl ) #/pl.sum()
        #cntr = pl.sum( 1./weightsr ) #/pl.sum()
        #cntc = pl.sum( 1./weightsc ) #/pl.sum()
        #print pixelsl, pixelsr, pixelsc
        #print cntl ,cntr, cntw #, weights, weightsl
        
        #areas = areaps*ps*ps*degtosterad  # area in sterad in which each source can be detected
        #areasum = pl.sum(1./areas)
        
        
        #areas =   # area in sterad in which each source can be detected
        areasum1 = pl.sum(1./(areafscales*Area_sttot))
        
        
        #areas =   # area in sterad in which each source can be detected
        areasum2 = pl.sum(1./(areafs*Area_sttot))
        
        #pl.figure(1)
        #try:
            #pl.hist(areaps, bins=20)
        #except:
            #print areaps
        pl.figure(2)
        pl.semilogx(fluxs, fluxs/peaks, 'o')
        mean_weight = pl.mean(1./areafs)
        
        nbin_err = cnt_err
        dnbin = bin_widths[i]
        counts_low[ i ] = ((cnt/Area_stl)/dnbin)*((bin_centres[i])**2.5)
        counts_high[ i ] = ((cnt/Area_str)/dnbin)*((bin_centres[i])**2.5)
        counts[ i ] = ((cnt/Area_stc)/dnbin)*((bin_centres[i])**2.5)
        counts_areascale[ i ] = ((cnt/Area_stc_areasclae)/dnbin)*((bin_centres[i])**2.5)
        #counts_weight[ i ] = ((cntw)/dnbin)*((bin_centres[i])**2.5)
        
        #counts_weight[ i ] = (areasum/dnbin)*((bin_centres[i])**2.5)
        counts_weight1[ i ] = (areasum1/dnbin)*((bin_centres[i])**2.5)
        counts_weight2[ i ] = (areasum2/dnbin)*((bin_centres[i])**2.5)
        
        #print cntl, cnt, cntr
        #print ((cntl/dnbin)/Area_stl)*((bin_centres[i])**2.5),((cnt/dnbin)/Area_stc)*((bin_centres[i])**2.5), ((cntr/dnbin)/Area_str)*((bin_centres[i])**2.5) 
        #print counts_low[ i ], counts[ i ], counts_high[ i ] 
        
        #counts_cor[i] = counts_weight1[i] #*corF[i]
        #count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
        #count_cor_errors[ i ] = pl.sqrt( (counts_cor[i]*count_errors[ i ]/counts[i])**2. +(counts_cor[i]*corFe[i]/corF[i])**2. )
        count_errors_l[ i ] = ((cnt_l_err/dnbin)/Area_stc)*((bin_centres[i])**2.5)
        count_errors_u[ i ] = ((cnt_u_err/dnbin)/Area_stc)*((bin_centres[i])**2.5)
        
        #count_errors_l[ i ] = pl.sqrt( (counts_cor[i]*count_errors_l[ i ]/counts[i])**2. +(counts_cor[i]*corFe[i]/corF[i])**2. )
        #count_errors_u[ i ] = pl.sqrt( (counts_cor[i]*count_errors_u[ i ]/counts[i])**2. +(counts_cor[i]*corFe[i]/corF[i])**2. )
        
        Areal = np.round(Areal, 2)
        Arear = np.round(Arear, 2)
        
        if Areal != Arear:
            print r'${b1:7.2f}-{b2:7.2f}$ & ${bcen:7.2f}$ & ${rcnt:3.0f}_{{-{rcntlerr:2.0f}}}^{{+{rcntuerr:2.0f}}}$ & ${areal:5.2f}-{arear:5.2f}$ & ${area:5.2f}$ & ${mw:.2f}$& ${corfdr:3.1f}$ & ${corc:3.1f}$  & ${corf:3.1f} \pm {corfe:3.1f}$ & ${cnt:4.0f}_{{-{cntlerr:3.0f}}}^{{+{cntuerr:4.0f}}}$ \\[+1.5pt]'.format(b1=bins[i]*scaleF, b2=bins[i+1]*scaleF, bcen=bin_centres[ i ]*scaleF, rcnt=cnt, rcntlerr=cnt_l_err, rcntuerr=cnt_u_err, areal=Areal, arear=Arear, area=Areac , mw=mean_weight , corfdr=counts_corfdr[i], corc=counts_corcomp[i],  corf=counts_corres[i], corfe=counts_corres_err[i], cnt=counts_corfdr[i]*counts_corcomp[i]*counts_corres[i]*counts_weight1[i], cntlerr=count_errors_l[i], cntuerr=count_errors_u[i])
        else:
            print r'${b1:7.2f}-{b2:7.2f}$ & ${bcen:7.2f}$ & ${rcnt:3.0f}_{{-{rcntlerr:2.0f}}}^{{+{rcntuerr:2.0f}}}$ &  $  \ldots $  & ${area:5.2f}$ & ${mw:.2f}$& ${corfdr:3.1f}$ & ${corc:3.1f}$ & ${corf:3.1f} \pm {corfe:3.1f}$ & ${cnt:4.0f}_{{-{cntlerr:3.0f}}}^{{+{cntuerr:4.0f}}}$ \\[+1.5pt]'.format(b1=bins[i]*scaleF, b2=bins[i+1]*scaleF, bcen=bin_centres[ i ]*scaleF, rcnt=cnt, rcntlerr=cnt_l_err, rcntuerr=cnt_u_err, areal=Areal, arear=Arear, area=Areac , mw=mean_weight, corfdr=counts_corfdr[i], corc=counts_corcomp[i] , corf=counts_corres[i], corfe=counts_corres_err[i], cnt=counts_corfdr[i]*counts_corcomp[i]*counts_corres[i]*counts_weight1[i], cntlerr=count_errors_l[i], cntuerr=count_errors_u[i])
        #if Areal != Arear:
            #print r'${b1:7.2f}-{b2:7.2f}$ & ${bcen:7.2f}$ & ${rcnt:3.0f}_{{-{rcntlerr:2.0f}}}^{{+{rcntuerr:2.0f}}}$ & ${areal:5.2f}-{arear:5.2f}$ & ${area:5.2f}$ & ${mw:.2f}$& ${corc:3.1f}$  & ${corf:3.1f} \pm {corfe:3.1f}$ & ${cnt:4.0f}_{{-{cntlerr:3.0f}}}^{{+{cntuerr:4.0f}}}$ \\[+1.5pt]'.format(b1=bins[i]*scaleF, b2=bins[i+1]*scaleF, bcen=bin_centres[ i ]*scaleF, rcnt=cnt, rcntlerr=cnt_l_err, rcntuerr=cnt_u_err, areal=Areal, arear=Arear, area=Areac , mw=mean_weight , corc=counts_corcomp[i],  corf=counts_corres[i], corfe=counts_corres_err[i], cnt=counts_corfdr[i]*counts_corcomp[i]*counts_corres[i]*counts_weight1[i], cntlerr=count_errors_l[i], cntuerr=count_errors_u[i])
        #else:
            #print r'${b1:7.2f}-{b2:7.2f}$ & ${bcen:7.2f}$ & ${rcnt:3.0f}_{{-{rcntlerr:2.0f}}}^{{+{rcntuerr:2.0f}}}$ &  $  \ldots $  & ${area:5.2f}$ & ${mw:.2f}$& ${corc:3.1f}$ & ${corf:3.1f} \pm {corfe:3.1f}$ & ${cnt:4.0f}_{{-{cntlerr:3.0f}}}^{{+{cntuerr:4.0f}}}$ \\[+1.5pt]'.format(b1=bins[i]*scaleF, b2=bins[i+1]*scaleF, bcen=bin_centres[ i ]*scaleF, rcnt=cnt, rcntlerr=cnt_l_err, rcntuerr=cnt_u_err, areal=Areal, arear=Arear, area=Areac , mw=mean_weight, corc=counts_corcomp[i] , corf=counts_corres[i], corfe=counts_corres_err[i], cnt=counts_corfdr[i]*counts_corcomp[i]*counts_corres[i]*counts_weight1[i], cntlerr=count_errors_l[i], cntuerr=count_errors_u[i])
    
        #sys.exit()
    
    
    
    fig, ax1 = pp.paper_double_ax()
    pl.minorticks_on()

    #ax1.xaxis.set_major_locator( pl.MaxNLocator(nbins=6, prune='lower') )
    #ax1.xaxis.set_minor_locator(MultipleLocator(0.25))
    #ax1.yaxis.set_major_locator( mpl.ticker.LogLocator(numticks=5))
    #ax1.yaxis.set_major_formatter( mpl.ticker.LogFormatter(base=10.0, labelOnlyBase=False))
    #ax1.xaxis.set_major_locator( mpl.ticker.LogLocator(numticks=5))
    #ax1.xaxis.set_major_formatter( mpl.ticker.LogFormatter(base=10.0, labelOnlyBase=False))
    #subs = [1.0, 2.0, 5.0]  # ticks to show per decade
    #ax1.xaxis.set_major_locator(mpl.ticker.LogLocator(subs=subs)) #set the ticks position
    #ax.xaxis.set_major_formatter(mpl.ticker.NullFormatter())   # remove the major ticks
    #ax1.xaxis.set_major_formatter(mpl.ticker.LogFormatter())  #add the custom ticks
    
    #ax1.loglog(bin_centres, counts_cor, 'bo')


    othercol = pl.cm.Set1(np.linspace(0, 1, 10))
    othersymb = ['x','*','o','^','v','>','<','d','p','*','h','H','1','2','3','4','D']
    othersymbi = 0

    Mauch_frq, Mauch_bins, Mauch_counts, Mauch_count_errors = get_Mauch_counts(nu, alpha=alpha)
    errMauch = ax1.errorbar(Mauch_bins*scaleF, Mauch_counts, Mauch_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='325 MHz, Mauch et al. 2013')
    othersymbi += 1
    

    
    gmrt_frq, gmrt_bins, gmrt_counts, gmrt_count_errors = get_gmrt_counts(nu, alpha=alpha, scale=False)
    errWilliams = ax1.errorbar(gmrt_bins*scaleF, gmrt_counts, gmrt_count_errors, fmt='k.', color='g', ms = 7, label='153 MHz, Williams et al. 2013')
    #ax1.errorbar(save_gmrt_bins*scaleF, save_gmrt_counts, save_gmrt_count_errs, fmt='ks', color='gray', ms = 7, label='Williams et al. 2013')
    counts=counts_cor
    
    print bin_centres
    print counts
    i1 = pl.sum(bin_centres<10.e-3)
    i2 = pl.sum(bin_centres<400.e-3)
    print i1,i2
    print pl.log10(bin_centres[i1:i2]), pl.log10(counts[i1:i2])
    try:
        pm, pc = np.polyfit(pl.log10(bin_centres[i1:i2]), pl.log10(counts[i1:i2]), 1)
        p0 = [pc,pm]
        print p0
    except:
        print "error polyfit"
        p0 = [0,1]
    xfit = pl.log10(bin_centres[i1:i2])
    yfit = pl.log10(counts[i1:i2])
    res = lambda p, xfit, yfit : (yfit-Tpoly(p, xfit))
    try:
        (p, cov, info, mesg, flag) = leastsq(res, p0, args=(xfit, yfit), full_output=True)
    except:
        cov = None
        p = p0
    if cov != None:
        chisq = sum(info["fvec"]*info["fvec"])
        dof = len(info["fvec"])-len(p)
        ep = np.array([np.sqrt(cov[i,i]*chisq/dof) for i in range(len(p))])
    else:
        p, ep = [np.nan, np.nan], [np.nan, np.nan]
    
    print 'power law intercept: %.3f \pm %.3f' %(p[0],ep[0])
    print 'power law slope: %.3f \pm %.3f' %(p[1],ep[1])
    
        
    #IC10_frq, IC10_bins, IC10_counts, IC10_count_errors = get_IC10_counts()
    #IC10_bins = IC10_bins*(nu/IC10_frq)**alpha
    #ax1.errorbar(IC10_bins, IC10_counts, IC10_count_errors, fmt='r^', ms = 5, mfc = 'r',label='Ishwara-Chandra et al. 2010')
    
    
    w_frq, w_bins, w_counts, w_count_errors = get_deZotti_counts(nu, alpha=alpha, value='all')
    errZotti = ax1.errorbar(w_bins*scaleF, w_counts, w_count_errors, color='gray', fmt='.', ms = 5,label='de Zotti et al. 2010')
    
    
    
    f1400, S1400, C1400, b1400  = get_Bondi_counts(nu, alpha=1)
    f, S0, C0,b0  = get_Bondi_counts(nu, alpha=-0.8)
    f, S1, C1,b1  = get_Bondi_counts(nu, alpha=-0.5)
    ## alpha range
    #C0 = 1.e0
    #S0 = 1.e1
    #a1 = -0.8
    #a2 = -0.5
    #dx = 1.5*(a1-a2)*np.log10(S0*nu/1400.)
    #dy = (a1-a2)*np.log10(S0*nu/1400.)
    #ax1.quiver(scaleF*S0[0], C0[0], scaleF*(S1[0]-S0[0]), C1[0]-C0[0], units='xy')
    #ax1.arrow(scaleF*S0[0], C0[0], scaleF*(S1[0]-S0[0]), C1[0]-C0[0], units='xy')
    #ax1.plot(scaleF*S0[0], C0[0], 'b', ms=20)
    #ax1.plot(scaleF*S1[0], C1[0], 'b', ms=20)
    #ax1.plot([scaleF*S0[0],scaleF*S1[0]], [C0[0],C1[0]], 'b', lw=3)
    #ax1.arrow(scaleF*S1400[0], C1400[0], scaleF*(S0[0]-S1400[0]), C0[0]-C1400[0], color='lime',lw=1)
    #ax1.arrow(scaleF*S1400[0], C1400[0], scaleF*(S1[0]-S1400[0]), C1[0]-C1400[0], color='magenta',lw=1)
    #ax1.arrow(scaleF*S1400[-1], C1400[-1], scaleF*(S0[-1]-S1400[-1]), C0[-1]-C1400[-1], color='lime',lw=1)
    #ax1.arrow(scaleF*S1400[-1], C1400[-1], scaleF*(S1[-1]-S1400[-1]), C1[-1]-C1400[-1], color='magenta',lw=1)
    
    ax1.plot([scaleF*S0[0], scaleF*S0[-1]], [C0[0], C0[-1]], color='blue',lw=2,alpha=0.5)
    ax1.plot([scaleF*S1[0], scaleF*S1[-1]], [C1[0], C1[-1]], color='blue',lw=2,alpha=0.5)
    ax1.plot([scaleF*S0[0], scaleF*S1[0]], [C0[0], C1[0]], color='blue', ls=':',lw=2)
    ax1.plot([scaleF*S0[-1], scaleF*S1[-1]], [C0[-1], C1[-1]], color='blue', ls=':',lw=2)#, length_includes_head=True, head_width=8, head_length=6)
    
    #ax1.arrow(scaleF*S0[0], C0[0], scaleF*(S1[0]-S0[0]), C1[0]-C0[0], color='gray',lw=1)#, length_includes_head=True, head_width=8, head_length=6)
    #ax1.arrow(scaleF*S0[-1], C0[-1], scaleF*(S1[-1]-S0[-1]), C1[-1]-C0[-1], color='gray',lw=1)#, length_includes_head=True, head_width=8, head_length=6)
    
    
    #ax1.errorbar(scaleF*S1400, C1400, b1400, color='gray', fmt='o', ms = 5)
    #ax1.errorbar(scaleF*S1, C1, b1, color='red', fmt='o', ms = 5)
    #ax1.errorbar(scaleF*S0, C0, b0, color='magenta', fmt='o', ms = 5)
    
    
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_agnonly() # AGN
    skadslineAGN151 = ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = '--',label='Wilman et al. 2008 - AGN')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=0, less=False) # SF
    skadslineSF151 = ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = ':',label='Wilman et al. 2008 - SF')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1,2],[1,2,3,4]]) # All
    skadslineAll151 = ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = '-',label='Wilman et al. 2008 - All')
    
    
    
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_agnonly(freq=610.) # AGN
    skadslineAGN, = ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = '--', label='Wilman et al. 2008 - AGN')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=0, less=False, freq=610.) # SF
    skadslineSF, = ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = ':', label='Wilman et al. 2008 - SF')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1,2],[1,2,3,4]], freq=610.) # All
    skadslineAll, = ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = '-', label='Wilman et al. 2008 - All')


    
    count_errors_u[count_errors_u>=counts_weight1] = 0.99*counts_weight1[count_errors_u>=counts_weight1]
    count_errors_l[count_errors_l>=counts_weight1] = 0.99*counts_weight1[count_errors_l>=counts_weight1]
    
    #ax1.errorbar(bin_centres*scaleF, counts_areascale, [count_errors_l,count_errors_u], [bin_erru*scaleF,bin_errl*scaleF], fmt='k.', ms = 10, label='This work peak area')
    
    #ax1.errorbar(bin_centres*scaleF, counts_low, fmt='kv', label='This work low')
    #ax1.errorbar(bin_centres*scaleF, counts, fmt='ko', label='This work mid')
    #ax1.errorbar(bin_centres*scaleF, counts_high, fmt='k^', label='This work high')
    #ax1.loglog(bin_centres*scaleF, counts, 'kx', label='This work cen')
    #ax1.loglog(bin_centres*scaleF, counts_area, 'ks', label='This work area')
    #ax1.loglog(bin_centres*scaleF, counts_weight, 'k+', label='This work weight')
    #ax1.errorbar(bin_centres*scaleF, counts_weight1, [count_errors_l,count_errors_u], [bin_errl*scaleF,bin_erru*scaleF], fmt='ko', markerfacecolor='none', label='This work')
    thisnocor, = ax1.plot(bin_centres*scaleF, counts_weight1,'ko', markerfacecolor='none', label='This work no correction')
    #count_errors_l_cor = np.max(np.array((count_errors_l, counts_corres*counts_weight1-counts_weight1)), axis=0)
    #count_errors_u_cor = np.max(np.array((count_errors_u, counts_weight1-counts_corres*counts_weight1)), axis=0)
    
    counts_cor_final = counts_corfdr*counts_corcomp*counts_corres*counts_weight1
    
    counts_err_correslow = counts_cor_final*np.sqrt((counts_corres_err/counts_corres)**2. + (count_errors_l/counts_weight1)**2.)
    counts_err_corresup = counts_cor_final*np.sqrt((counts_corres_err/counts_corres)**2. + (count_errors_u/counts_weight1)**2.)
    
    thisdata = ax1.errorbar(bin_centres*scaleF, counts_cor_final, [counts_err_correslow,counts_err_corresup], [bin_errl*scaleF,bin_erru*scaleF], markersize=10, fmt='ko', label='This work')
    #ax1.errorbar(bin_centres*scaleF, counts_corres*counts_weight1, [count_errors_l_cor,count_errors_u_cor], [bin_errl*scaleF,bin_erru*scaleF], fmt='ko', label='This work')
    #ax1.plot(bin_centres*scaleF, counts_correslow*counts_weight1, 'kv', label='This work weight')
    #ax1.plot(bin_centres*scaleF, counts_correshigh*counts_weight1, 'k^', label='This work weight')
    #ax1.loglog(bin_centres*scaleF, counts_weight2, 'm+', label='This work reweight scale peak')
    #ax1.loglog(bin_centres*scaleF, 2*counts_weight1, 'y+', label='This work weight x2')
    #ax1.loglog(bin_centres, counts_low, 'gv')
    #ax1.loglog(bin_centres, counts_high, 'c^')
    #ax1.loglog(bin_centres, counts_weight, 'ro')
    #ax1.errorbar(bin_centres, counts_low, count_errors, fmt='ko', ms = 7, label='This work')
    #ax1.errorbar(bin_centres, counts_cor, count_cor_errors, fmt='ko', ms = 7, label='This work cor')
    
    
    src_counts_nocor = [bin_centres*scaleF, counts_weight1, count_errors_l,count_errors_u]
    np.save('%s_source_counts_nocor.npy' %(name), src_counts_nocor)
    src_counts = [bin_centres*scaleF, counts_cor_final, counts_err_correslow,counts_err_corresup]
    np.save('%s_source_counts.npy' %(name), src_counts)
    
    src_counts=np.asarray(src_counts)
    src_counts = src_counts.T
    for i in range(len(src_counts)):
        print '%f\t%f\t%f\t%f' %(src_counts[i,0],src_counts[i,1], src_counts[i,2], src_counts[i,3])
        
    src_counts_nocor=np.asarray(src_counts_nocor)
    src_counts_nocor = src_counts_nocor.T
    for i in range(len(src_counts_nocor)):
        print '%f\t%f\t%f\t%f' %(src_counts_nocor[i,0],src_counts_nocor[i,1], src_counts_nocor[i,2], src_counts_nocor[i,3])
    
    dummy = ax1.plot()
    
    y1,y2 = ax1.get_ylim()
    
    y1 ,y2 = 0.15,9000
    #ax1.vlines(PSNRLIM_DET*rms10*scaleF, y1,y2, linestyles='dashed')
    #ax1.vlines(PSNRLIM_DET*rms50*scaleF, y1,y2, linestyles='dashed')
    #ax1.vlines(PSNRLIM_DET*rms90*scaleF, y1,y2, linestyles='dashed')
    
    #pl.title(' Source Counts ')
    pl.ylabel(r'$S^{5/2} dN / dS$ [Jy$^{3/2}$ sr$^{-1}$]')
    pl.xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='small'
    #pl.legend(loc=4,numpoints=1,frameon=False,ncol=2)
    # Create a legend for the first line.
    first_legend = pl.legend([skadslineAll, skadslineAGN, skadslineSF, errMauch, errWilliams, thisnocor, thisdata] ,['Wilman et al. 2008 - All', 'Wilman et al. 2008 - AGN', 'Wilman et al. 2008 - SF', '325 MHz, Mauch et al. 2013', '153 MHz, Williams et al. 2013','This work no correction', 'This work'], loc='lower right', numpoints=1, fontsize='small', frameon = True)
    frame = first_legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')
    #first_legend = pl.legend([skadslineAll, skadslineAGN, skadslineSF, errMauch, errWilliams, thisnocor, thisdata] ,['Wilman et al. 2008 - All', 'Wilman et al. 2008 - AGN', 'Wilman et al. 2008 - SF', '325 MHz, Mauch et al. 2013', '153 MHz, Williams et al. 2013','This work no correction', 'This work'], loc='lower right', frameon=False, numpoints=1, fontsize='small')
    # Add the legend manually to the current Axes.
    ax = pl.gca().add_artist(first_legend)
    # Create another legend for the second line.
    #second_legend = pl.legend([errWhite, errCiliegi, errRichards, errPrandoni, errHopkins, errBondi, errSeymour, errSimpson], ['White et al. 1997', 'Ciliegi et al. 1999', 'Richards 2000', 'Prandoni et al. 2001', 'Hopkins et al. 2002', 'Bondi et al. 2003', 'Seymour et al. 2004', 'Simpson et al. 2006'], loc='upper left', frameon=True, numpoints=1, title='1.4 GHz Scaled Counts', fontsize='small')
    second_legend = pl.legend([errZotti], ['de Zotti et al. 2010'], loc='upper left', frameon=True, numpoints=1, title='1.4 GHz Scaled Counts', fontsize='small')
    frame = second_legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')
    
    
    
    
    #pl.minorticks_on()
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))
    
    ax1.set_xlim(0.15,9000)
    ax1.set_ylim(2, 9000)
    
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    
    pp.fig_save_many(fig, '%s_source_counts' %(name))
    

######

    fig, ax1 = pp.paper_double_ax()
    pl.minorticks_on()

    #othercol = pl.cm.Set1(np.linspace(0, 1, 10))
    #othersymb = ['x','*','o','^','v','>','<','d','p','*','h','H','1','2','3','4','D']
    #othersymbi = 0

    #Mauch_frq, Mauch_bins, Mauch_counts, Mauch_count_errors = get_Mauch_counts(nu, alpha=alpha)
    #errMauch = ax1.errorbar(Mauch_bins*scaleF, Mauch_counts, Mauch_count_errors, color=othercol[othersymbi], fmt=othersymb[othersymbi], ms = 5,label='325 MHz, Mauch et al. 2013')
    #othersymbi += 1
    
    #gmrt_frq, gmrt_bins, gmrt_counts, gmrt_count_errors = get_gmrt_counts(nu, alpha=alpha, scale=False)
    #errWilliams = ax1.errorbar(gmrt_bins*scaleF, gmrt_counts, gmrt_count_errors, fmt='k.', color='g', ms = 7, label='153 MHz, Williams et al. 2013')
    #counts=counts_cor
    
    #print bin_centres
    #print counts
    #i1 = pl.sum(bin_centres<10.e-3)
    #i2 = pl.sum(bin_centres<400.e-3)
    #print i1,i2
    #print pl.log10(bin_centres[i1:i2]), pl.log10(counts[i1:i2])
    #try:
        #pm, pc = np.polyfit(pl.log10(bin_centres[i1:i2]), pl.log10(counts[i1:i2]), 1)
        #p0 = [pc,pm]
        #print p0
    #except:
        #print "error polyfit"
        #p0 = [0,1]
    #xfit = pl.log10(bin_centres[i1:i2])
    #yfit = pl.log10(counts[i1:i2])
    #res = lambda p, xfit, yfit : (yfit-Tpoly(p, xfit))
    #try:
        #(p, cov, info, mesg, flag) = leastsq(res, p0, args=(xfit, yfit), full_output=True)
    #except:
        #cov = None
        #p = p0
    #if cov != None:
        #chisq = sum(info["fvec"]*info["fvec"])
        #dof = len(info["fvec"])-len(p)
        #ep = np.array([np.sqrt(cov[i,i]*chisq/dof) for i in range(len(p))])
    #else:
        #p, ep = [np.nan, np.nan], [np.nan, np.nan]
    
    #print 'power law intercept: %.3f \pm %.3f' %(p[0],ep[0])
    #print 'power law slope: %.3f \pm %.3f' %(p[1],ep[1])
    
    #w_frq, w_bins, w_counts, w_count_errors = get_deZotti_counts(nu, alpha=alpha, value='all')
    #errZotti = ax1.errorbar(w_bins*scaleF, w_counts, w_count_errors, color='gray', fmt='.', ms = 5,label='de Zotti et al. 2010')
    
    #f1400, S1400, C1400, b1400  = get_Bondi_counts(nu, alpha=1)
    #f, S0, C0,b0  = get_Bondi_counts(nu, alpha=-0.8)
    #f, S1, C1,b1  = get_Bondi_counts(nu, alpha=-0.5)

    #ax1.plot([scaleF*S0[0], scaleF*S0[-1]], [C0[0], C0[-1]], color='blue',lw=2,alpha=0.5)
    #ax1.plot([scaleF*S1[0], scaleF*S1[-1]], [C1[0], C1[-1]], color='blue',lw=2,alpha=0.5)
    #ax1.plot([scaleF*S0[0], scaleF*S1[0]], [C0[0], C1[0]], color='blue', ls=':',lw=2)
    #ax1.plot([scaleF*S0[-1], scaleF*S1[-1]], [C0[-1], C1[-1]], color='blue', ls=':',lw=2)#, length_includes_head=True, head_width=8, head_length=6)
   
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_agnonly() # AGN
    skadslineAGN151 = ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = '--',label='Wilman et al. 2008 - AGN')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=0, less=False) # SF
    skadslineSF151 = ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = ':',label='Wilman et al. 2008 - SF')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1,2],[1,2,3,4]]) # All
    skadslineAll151 = ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = '-',label='Wilman et al. 2008 - All')
    
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_agnonly(freq=610.) # AGN
    skadslineAGN, = ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = '--', label='Wilman et al. 2008 - AGN')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=0, less=False, freq=610.) # SF
    skadslineSF, = ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = ':', label='Wilman et al. 2008 - SF')
    skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1,2],[1,2,3,4]], freq=610.) # All
    skadslineAll, = ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = '-', label='Wilman et al. 2008 - All')

    count_errors_u[count_errors_u>=counts_weight1] = 0.99*counts_weight1[count_errors_u>=counts_weight1]
    count_errors_l[count_errors_l>=counts_weight1] = 0.99*counts_weight1[count_errors_l>=counts_weight1]
    #thisnocor, = ax1.plot(bin_centres*scaleF, counts_weight1,'ko', markerfacecolor='none', label='This work no correction')
    
    counts_err_correslow = counts_weight1*np.sqrt((counts_corres_err/counts_corres)**2. + (count_errors_l/counts_weight1)**2.)
    counts_err_corresup = counts_weight1*np.sqrt((counts_corres_err/counts_corres)**2. + (count_errors_u/counts_weight1)**2.)
    
    thisdata = ax1.errorbar(bin_centres*scaleF, counts_corres*counts_weight1, [counts_err_correslow,counts_err_corresup], [bin_errl*scaleF,bin_erru*scaleF], markersize=10, fmt='ko', label='This work')
    
    dummy = ax1.plot()
    
    
    y1 ,y2 = 0.15,9000
    
    pl.ylabel(r'$S^{5/2} dN / dS$ [Jy$^{3/2}$ sr$^{-1}$]')
    pl.xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='small'
    #pl.legend(loc=4,numpoints=1,frameon=False,ncol=2)
    # Create a legend for the first line.
    first_legend = pl.legend([skadslineAll, skadslineAGN, skadslineSF, thisdata] ,['Wilman et al. 2008 - All', 'Wilman et al. 2008 - AGN', 'Wilman et al. 2008 - SF', 'This work'], loc='lower right', numpoints=1, fontsize='small', frameon = True)
    frame = first_legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')
    #first_legend = pl.legend([skadslineAll, skadslineAGN, skadslineSF, errMauch, errWilliams, thisnocor, thisdata] ,['Wilman et al. 2008 - All', 'Wilman et al. 2008 - AGN', 'Wilman et al. 2008 - SF', '325 MHz, Mauch et al. 2013', '153 MHz, Williams et al. 2013','This work no correction', 'This work'], loc='lower right', frameon=False, numpoints=1, fontsize='small')
    # Add the legend manually to the current Axes.
    ax = pl.gca().add_artist(first_legend)
    # Create another legend for the second line.
    #second_legend = pl.legend([errWhite, errCiliegi, errRichards, errPrandoni, errHopkins, errBondi, errSeymour, errSimpson], ['White et al. 1997', 'Ciliegi et al. 1999', 'Richards 2000', 'Prandoni et al. 2001', 'Hopkins et al. 2002', 'Bondi et al. 2003', 'Seymour et al. 2004', 'Simpson et al. 2006'], loc='upper left', frameon=True, numpoints=1, title='1.4 GHz Scaled Counts', fontsize='small')
    #second_legend = pl.legend([errZotti], ['de Zotti et al. 2010'], loc='upper left', frameon=True, numpoints=1, title='1.4 GHz Scaled Counts', fontsize='small')
    #frame = second_legend.get_frame()
    #frame.set_facecolor('white')
    #frame.set_edgecolor('white')
    
    ax1.set_xlim(0.15,9000)
    ax1.set_ylim(2, 9000)
    
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    
    pp.fig_save_many(fig, '%s_source_counts_simple' %(name))
    

######



    
    
    #fig = pl.figure(figsize=(7.0,5.0))
    #ax1 = pl.subplot(1, 1, 1)
    #ax1.loglog(bin_centres, counts, 'k')
    #ax1.errorbar(bin_centres, counts, count_errors, fmt='ko', ms = 7)
    
    ##ax1.set_xlim(2e-3, 3.)
    ##ax1.set_ylim(1.e1, 1.e4)
    ##pl.title(' Source Counts ')
    #pl.ylabel(r'$S^{5/2} dN / dS$ [Jy$^{3/2}$ sr$^{-1}$]')
    #pl.xlabel(r'%i MHz flux [Jy]' %(nu))
    
    #pl.minorticks_on()
    #pl.savefig('%s_source_counts_specind.png' %(name))
    #pl.savefig('%s_source_counts_specind.eps' %(name))


    #src_counts = [bin_centres, counts_cor, count_errors_l,count_errors_u]
    #np.save('%s_source_counts.npy' %(name), src_counts)
    savefile = open('%s_source_counts.dat' %(name),'w')
    savefile.write('#bin_centre\tcorrected_counts\tcount_err_l\tcount_err_u\n')
    for i in range(nbins):
        s = '%f\t%f\t%f\t%f\n' %(bin_centres[i], counts_cor[i], count_errors_l[i],count_errors_u[i])
    savefile.write(s)
    savefile.close()
    
    plot_src_counts_all()

    #data = pl.array(pl.loadtxt('skads_151flux_1mJy_type_z.result',delimiter=',',skiprows=1)).transpose()
    #return



#def main():
    
    #parser = argparse.ArgumentParser() 
    #parser.add_argument('fits_cat', help="Name of fits srouce catalogue")
    #parser.add_argument('rmsfits', help="Name of rms fits image for area")
    #parser.add_argument('name', help="Identifying Name")
    #parser.add_argument('nu', help="Frequency")
    #parser.add_argument('--unit', default='mJy', help="(default mJy)")
    #parser.add_argument('-c','--clobber', action='store_true', default=False, help="(default false)")
    
    #args = parser.parse_args()
    
    #fits_cat = args.fits_cat
    #rmsfits = args.rmsfits
    #name = args.name
    #nu = float(args.nu)
    #plot_source_counts(fits_cat, rmsfits, name, nu, corFlux=1., unit=args.unit, clobber_areas=args.clobber)



#if __name__ == "__main__":
    #main()




'''
figure()
n,b,p = hist(flux[flux>2e-3], range=(2e-3,1), bins=100 )
bw = b[1:] - b[:-1]
dnds = n/bw

figure()
loglog(b[:-1], dnds)

x = b[:-1]
y = dnds

logx = np.log10(x)
logy = np.log10(y)
logx = logx[np.isfinite(logy)]
logy = logy[np.isfinite(logy)]

p = np.polyfit(logx,logy,1)

A = 10**(p[1])
beta = p[0]

py = A*x**beta



plot(x,py)


b1 = 0.000607
b2 = 1.258
A = 51.18
beta = -1.436
Nperbin = 50.
bins=[b1]
while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))

figure()
n1,b1,p1 = hist(flux, bins=bins )
'''
