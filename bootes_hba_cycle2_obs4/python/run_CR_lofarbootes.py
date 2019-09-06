import matplotlib as mpl
mpl.use('Agg')
#mpl.rc_file('/home/wwilliams/.config/matplotlib/matplotlibrc')  # <-- the file containing your settings
#mpl.rc_file('/local/phd/gmrt_bootes/a2256_mosaic/fits/matplotlibrc') 
from fits_util import load_fits
from CR_core_lofar import *
import sys

CLOBBER=True
CLOBBER=False

NsimP = 1
NsimE = 1


if len(sys.argv) > 1:
    simi = int(sys.argv[1])
    #PSIMS = range(NsimP)
    #ESIMS = range(NsimE)
    NsimP = 7
    PSIMS = range(0,NsimP)
    #PSIMS = [simi]
    ESIMS = [simi]
    
else:
    NsimP = 7
    NsimE = 7
    PSIMS = range(0,NsimP)
    ESIMS = range(0,NsimE)
    

#FIELDNAME = 'BOOTES_LOFAR'
#Nsim = 25


#find_sources_pybdsm.py --incl_empty  --ungroup --group_tol 10.0 --rmsbox 160,50 --rmsboxbright 60,15 --detectimage mosaic_v4t.pbcut.fits mosaic_v4t.pbcor.fits

find_sources_cmd = "python /home/wwilliams/scripts/find_sources_pybdsm.py  --incl_empty  --ungroup --group_tol 10.0 --rmsbox 160,50 --rmsboxbright 60,15 -d DETECT IMAGE"
#find_sources_cmd = "python ~/para/scripts/find_sources_pybdsm.py --atrous_do --ungroup --group_tol 10.0 --rmsbox 160,50 --rmsboxbright 80,20 -d DETECT IMAGE"

#PATH = '/data/lofar/wwilliams/bootes/bootes_hba_cycle2_obs4/mosaic_v4/CRtemp'




#SOURCElist = PATH + '../catalogue_v0.1/mosaic_v0.1_cutout.pbcor.fits.pybdsm.fits'
#RESIDUALmap = PATH + '../catalogue_v0.1/mosaic_v0.1_cutout.pbcor.fits.pybdsm.res.fits'
#RMSmap = PATH + '../catalogue_v0.1/mosaic_v0.1_cutout.pbcor.fits.pybdsm.rms.fits'
#INVERTmap = PATH + 'mosaic_v0.1_cutout.pbcor.inv.fits'

#facet = 's15'
#PATH = '/data1/wwilliams/hba_bootes_lc015/BOOTES_24_V4/'
#SIMPATH = PATH+'CR/{s}/'.format(s=facet)
#INmapd = PATH + '../mosaic_v2_split/{s}.pbcut.fits'.format(s=facet)
#INmap = PATH + '../mosaic_v2_split/{s}.pbcor.fits'.format(s=facet)
#SOURCElist = PATH + '../mosaic_v2_split/{s}.pbcor.fits.pybdsm.fits'.format(s=facet)
#PBmap = PATH + '../mosaic_v2_split/{s}.pb.fits'.format(s=facet)
#RESIDUALmap = PATH + '../mosaic_v2_split/{s}.pbcor.fits.pybdsm.res.fits'.format(s=facet)
#RMSmap = PATH + '../mosaic_v2_split/{s}.pbcor.fits.pybdsm.rms.fits'.format(s=facet)
#INVERTmap = PATH + '../mosaic_v2_split/{s}.pbcor.inv.fits'.format(s=facet)


#PATH = '/data1/wwilliams/hba_bootes_lc015/BOOTES_24_V4/pbcor_mosaic/v1/'
#SIMPATH = 'CR/'
#facet = 'mosaic_v1_cutout'
#INmapd =  '{s}.pbcut.fits'.format(s=facet)
#INmap =  '{s}.pbcor.fits'.format(s=facet)
#SOURCElist =  SIMPATH+'{s}.pbcor.fits.pybdsm.fits'.format(s=facet)
#PBmap =  SIMPATH+'{s}.pb.fits'.format(s=facet)
#RESIDUALmap =  SIMPATH+'{s}.pbcor.fits.pybdsm.res.fits'.format(s=facet)
#RMSmap =  SIMPATH+'{s}.pbcor.fits.pybdsm.rms.fits'.format(s=facet)
#INVERTmap =  '{s}.pbcor.inv.fits'.format(s=facet)

if 'uhppc58' in os.uname()[1]:
    PATH = '/local/wwilliams/projects/radio_imaging/bootes_hba_cycle2_obs4/mosaic_v4/'
    local = True
else:
    PATH = '/car-data/wwilliams/bootes/bootes_hba_cycle2_obs4/mosaic_v4/'
    PATH = '/data/lofar/wwilliams/bootes/bootes_hba_cycle2_obs4/mosaic_v4/'
    local = False
SIMPATH = PATH+'CR/'
SIMPATH = PATH+'CRv2/'
SIMPATH = PATH+'CRv3/'
INmapd =   PATH+'mosaic_v4.pbcut.fits'
INmap =   PATH+'mosaic_v4.pbcor.fits'
PBmap =   PATH+'mosaic_v4.pb.fits'
SOURCElist = PATH+'../catalogue_v4/mosaic_v4.pbcor.pybdsm.corpos.corflux.fits'
RESIDUALmap = PATH+'mosaic_v4.pbcor.pybdsm.res.fits'
RMSmap =   PATH+'mosaic_v4.pbcor.pybdsm.rms.fits'
INVERTmap =  'mosaic_v4.pbcor.invert.fits'

INVERT_DIR = SIMPATH + 'invert/'
MODEL_DIR = SIMPATH + 'model/'

if not os.path.exists(SIMPATH):
    os.system("mkdir -p {p}".format(p=SIMPATH))

#os.chdir(SIMPATH)

# fix header info #
#fix_beam_in_header(INmap)
if (not os.path.exists(SOURCElist)): # or CLOBBER:
    #find_sources(PATH, INmap, SOURCElist, clobber=True, thresh_pix=5.0, thresh_isl=3.0, rms_box=(175,43), output_files=True)
    pwd = os.getcwd()
    os.chdir(SIMPATH)
    
    find_sources_cmd_map = find_sources_cmd.replace('DETECT', INmapd).replace('IMAGE',INmap)
    print find_sources_cmd_map
    os.system(find_sources_cmd_map)
  
    os.chdir( pwd )
    
#fix_beam_in_header(RESIDUALmap)



# get average rms #
sigma = rms_stat(RMSmap)

    
## assuming point sources ## Peak fluxes are 'cleaner' ##
t=load_fits(SOURCElist)
FLUXMAX = np.nanmax(t.Peak_flux)
FLUXMIN = np.nanmin(t.Peak_flux)
Fracmin = 0.25
Nobs = t.shape[0]
slope = get_flux_slope(t.Peak_flux,sigma)
print 'slope of peak flux distribution : %.3f' %(slope)


print 'sources in cat : %i' %(Nobs)
print 'flux range : %.3f - %.3f mJy' %(FLUXMIN*1000., FLUXMAX*1000.)

NSOURCES = Nobs*( Fracmin**slope*FLUXMIN**slope - FLUXMAX**slope )/( FLUXMIN**slope - FLUXMAX**slope )  
FLUXMIN = Fracmin*FLUXMIN

print 'sources to simulate : %f' %(NSOURCES)
print 'flux range : %.3f - %.3f mJy' %(FLUXMIN*1000., FLUXMAX*1000.)




Npoints = 15.
i1 = np.log10(0.0005)
i2 = np.log10(0.05)
irange = np.logspace(i1,i2,Npoints)
#Reliability(PBmap, INmapd, sigma, INVERT_DIR, INVERTmap, SOURCElist, irange, find_sources_cmd, savename=SIMPATH+'LOFAR_BOOTES', title='' )

#sys.exit()

## POINT SRC SIM ##
point_out_list = []
for simi in PSIMS:
    Sim = SimMap(INmap, SOURCElist, RESIDUALmap, PBmap, NSOURCES, FLUXMIN, FLUXMAX, MODEL_DIR, simi, alpha=slope, detect_thresh=5., isl_thresh=3.,  rms_box=(150,40), assocR=None, rootname='modelp', find_sources_args=find_sources_cmd)
    if not local:
        Sim.make_model_map(clobber=CLOBBER)
        Sim.find_sources_map(find_sources_cmd, clobber=CLOBBER)

    point_out_list.append(Sim.OutCatName)


p1 = 0.7   # 5
p2 = 3.0   # 1000
SNRrange = np.logspace(p1,p2,21)
PointFrac, PSNRenv = get_peak_total_env_frac(point_out_list, SOURCElist, SNRrange, savename=SIMPATH+'LOFAR_BOOTES', title='')

#sys.exit()
#CLOBBER = True

## EXT SRC SIM ##
matched_in_out_list = []
ext_out_list = []
for simi in ESIMS:
    Sim = SimMap(INmap, SOURCElist, RESIDUALmap, PBmap, NSOURCES, FLUXMIN, FLUXMAX, MODEL_DIR, simi, detect_thresh=5., isl_thresh=3.,  rms_box=(150,40), assocR=None, rootname='model', Fracpoint=PointFrac, find_sources_args=find_sources_cmd)
    if not local:
        Sim.make_model_map(clobber=CLOBBER)
        Sim.find_sources_map(find_sources_cmd, clobber=CLOBBER)
    Sim.match_catalogues(clobber=CLOBBER)

    matched_in_out_list.append(Sim.MatchInOutName)
    ext_out_list.append(Sim.OutCatName)

get_size_dist(ext_out_list, SOURCElist, SNRrange, savename=SIMPATH+'LOFAR_BOOTES', title='')

Npoints = 25.
i1 = np.log10(0.0005)
i2 = np.log10(0.05)
irange = np.logspace(i1,i2,Npoints)


diff_area_i = get_differential_area(RMSmap, irange)

Completeness(matched_in_out_list, sigma, irange, diff_area_i, savename=SIMPATH+'LOFAR_BOOTES',  title='', rmsfits=RMSmap)
