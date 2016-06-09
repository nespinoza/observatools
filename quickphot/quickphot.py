# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pyfits
import glob
import os
from photutils import CircularAperture,CircularAnnulus,aperture_photometry

EGAIN = 1.0 # electrons/ADU
ERON = 3.0 # electrons

def getApertureFluxes(subimg,x_cen,y_cen,Radius,sky_sigma,GAIN):
    apertures = CircularAperture([(x_cen,y_cen)],r=Radius)
    rawflux_table = aperture_photometry(subimg, apertures, \
            error=sky_sigma, effective_gain = GAIN)
    return rawflux_table['aperture_sum'][0],rawflux_table['aperture_sum_err'][0]

from scipy.ndimage.filters import gaussian_filter
def get_centroid(original_subimg):
    subimg = gaussian_filter(original_subimg,10)
    total = np.sum(subimg)
    x0 = np.sum(np.sum(subimg,axis=0)*np.arange(subimg.shape[0]))/np.sum(np.sum(subimg,axis=0))
    y0 = np.sum(np.sum(subimg,axis=1)*np.arange(subimg.shape[1]))/np.sum(np.sum(subimg,axis=1))
    return x0,y0

def read_dfile(fname):
    ff = open(fname,'r')
    im = []
    t = []
    f = []
    ferr = []
    while True:
        line = ff.readline()
        if line != '':
            d1,d2,d3,d4 = line.split()
            im.append(d1)
            t.append(np.double(d2))
            f.append(np.double(d3))
            ferr.append(np.double(d4))
        else:
            break
    ff.close()
    return im,t,f,ferr

def bin_data(t,rf_mag,rf_mag_err,n_bin = 10):
    times_bins = []
    mags_bins = []
    errors_bins = []
    for i in range(0,len(t),n_bin):
        times_bins.append(np.median(t[i:i+n_bin-1]))
        mags_bins.append(np.median(rf_mag[i:i+n_bin-1]))
        errors_bins.append(np.sqrt(np.sum(rf_mag_err[i:i+n_bin-1]**2))/np.double(n_bin))
    return times_bins,mags_bins,errors_bins

#############################################
filename = 'k11_quickphot.dat'
images = glob.glob('*.fits')
target_center = [2840,1520]
comp_center = [1870,1733]#[1080,660]
first_time = True
hs = 50
ap_rad = 25
#############################################

if not first_time:
   all_images,times,rel_fluxes,rel_fluxes_err = read_dfile(filename)
else:
    all_images = []
    times = []
    rel_fluxes = []
    rel_fluxes_err = []

if not os.path.exists('phot_images_target'):
   os.mkdir('phot_images_target')

if not os.path.exists('phot_images_comp'):
   os.mkdir('phot_images_comp')

for im in images:
    d,h = pyfits.getdata(im,header=True)
    if (h['FILTER'] == 'u') and (im not in all_images):
        all_images.append(im)
        time_ut = h['UT-TIME']
        hh,mm,ss = time_ut.split(':')
        hh = np.double(hh)
        mm = np.double(mm)
        ss = np.double(ss)
        if hh>0:
            time = hh + (mm/60.) + (ss/3600.)
        subimg = d[target_center[0]-hs:target_center[0]+hs,target_center[1]-hs:target_center[1]+hs]
        sky = np.median(subimg)
        subimg = subimg - sky
        imm = plt.imshow(subimg)
        imm.set_clim(0,1000)
        x0,y0 = get_centroid(subimg)
        plt.plot(x0,y0,'wx',markersize=15,alpha=0.5)
        circle = plt.Circle((x0,y0),ap_rad,color='white',fill=False)
        plt.gca().add_artist(circle)
        plt.savefig('phot_images_target/'+im+'.png')
        plt.close()
        f,ferr = getApertureFluxes(subimg,x0,y0,ap_rad,np.sqrt(sky),EGAIN)
        subimg = d[comp_center[0]-hs:comp_center[0]+hs,comp_center[1]-hs:comp_center[1]+hs]
        sky = np.median(subimg)
        subimg = subimg - sky
        imm = plt.imshow(subimg)
        imm.set_clim(0,1000)
        x0,y0 = get_centroid(subimg)
        plt.plot(x0,y0,'wx',markersize=15,alpha=0.5)
        circle = plt.Circle((x0,y0),ap_rad,color='white',fill=False)
        plt.gca().add_artist(circle)
        plt.savefig('phot_images_comp/'+im+'.png')
        plt.close()
        fc,ferrc = getApertureFluxes(subimg,x0,y0,ap_rad,np.sqrt(sky),EGAIN)
        times.append(time)
        rel_fluxes.append(f/fc)
        err = np.sqrt((ferr/fc)**2 + (f*ferrc/(fc**2))**2)
        rel_fluxes_err.append(err)      
plt.close()
print 'Writing...'
ff = open(filename,'w')
times = np.array(times)
rel_fluxes = np.array(rel_fluxes)
rel_fluxes_err = np.array(rel_fluxes_err)
rel_fluxes = rel_fluxes / np.median(rel_fluxes)
rel_fluxes_err = rel_fluxes_err / np.median(rel_fluxes)
for i in range(len(all_images)):
    ff.write(all_images[i]+'  '+str(times[i])+'  '+str(rel_fluxes[i])+'  '+str(rel_fluxes_err[i])+'\n')
ff.close()
plt.style.use('ggplot')
diff_mags = -2.5*np.log10(rel_fluxes)
err = (1.08574/rel_fluxes)*rel_fluxes_err
plt.errorbar(times,(1.-rel_fluxes)*1e6,yerr=rel_fluxes_err*1e6,fmt='o')
plt.xlabel('Time since mid-night (UT)')
plt.ylabel('Relative flux (ppm)')
plt.show()
idx = np.argsort(times)
times = times[idx][4:]
diff_mags = diff_mags[idx][4:]
err = err[idx][4:]
plt.errorbar(times,diff_mags-np.median(diff_mags[:20]),yerr=err,fmt='o',label = 'Data',alpha = 0.25)
tb,mb,mb_err = bin_data(times,diff_mags-np.median(diff_mags[:20]),err,n_bin=5)
plt.errorbar(tb,mb,yerr=mb_err,fmt='o',label = 'Binned data')
plt.legend()
plt.xlabel('Time since mid-night (UT)')
plt.ylabel('Differential magnitude')
plt.gca().invert_yaxis()
plt.show()
