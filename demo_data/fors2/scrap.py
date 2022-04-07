import astropy.io.fits as fits
import astroscrappy as asc
import copy
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np

files = glob.glob('tg*.fits')

n_iter   = 4
sig_clip = 9.
fs_mode  = 'convolve'
psfmod   = 'gaussx'
psfwidth = 3.9
obj_lim  = 2.

for f in files:
	with fits.open(f, mode='update') as han:
		head = han[0].header
		head['CRREMOVE'] = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S'), 'UT cosmic ray removal time'
		head['COMMENT'] = '  astroscrappy v1.0.8 used to remove cosmic rays'
		scihead = han[1].header
		ogdata = copy.deepcopy(han['SCI'].data)
		crmask, data = asc.detect_cosmics(copy.deepcopy(han['SCI'].data), 
		                                  gain=head['GAINMULT'], 
		                                  readnoise=head['RDNOISE'], 
		                                  niter=n_iter, 
		                                  sigclip=sig_clip,
		                                  fsmode=fs_mode,
		                                  psfmodel=psfmod,
		                                  psffwhm=3.9,
		                                  objlim=obj_lim)
		han['SCI'].data = data
		spec2Dhdu = fits.ImageHDU(ogdata)
		spec2Dhdu.header['EXTNAME'] = 'CRORIG_2D_SPEC'
		han['DQ'].data = (crmask*1) + han['DQ'].data
		han['DQ'].data[han['DQ'].data>1] = 1
		crmaskhdu = fits.ImageHDU(crmask*1)
		crmaskhdu.header['EXTNAME']  = '2D_CR_MASK'
		crmaskhdu.header['COMMENT']  = 'astroscrappy v1.0.8 used to remove cosmic rays'
		crmaskhdu.header['GAIN']     = head['GAINMULT']
		crmaskhdu.header['RDNOISE']  = head['RDNOISE']
		crmaskhdu.header['NITER']    = n_iter
		crmaskhdu.header['SIGCLIP']  = sig_clip
		crmaskhdu.header['FSMODE']   = fs_mode
		crmaskhdu.header['PSFMODEL'] = psfmod
		crmaskhdu.header['PSFFWHM']  = psfwidth
		crmaskhdu.header['OBJLIM']  = obj_lim
		han.append(spec2Dhdu)
		han.append(crmaskhdu)
		han.flush()
