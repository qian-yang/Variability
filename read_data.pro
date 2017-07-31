function mag2flux, mag
	flux = (10.0^((-0.4) * mag)) * (3.631*(1e-20)) ; erg s-1 Hz-1 cm-2
  return, flux
end

function magerr2fluxerr, mag, magerr, flux
	fluxerr = flux * 0.4 * alog(10.d) * magerr
	return, fluxerr
end

pro read_data, file, zydata = zydata
readcol, file, format = 'd, d, d', mjd, mag, magerr
ind = where((mjd gt 0) and (mag gt 0) and (magerr gt 0), nlc)
mjd = mjd[ind]
mag = transpose(mag[ind])
magerr = magerr[ind]
flux = mag2flux(mag)
fluxerr = magerr2fluxerr(mag, magerr, flux)
fluxivar = 1.0/(fluxerr^2)
cont_mean = mean(mag)
cont_mean_err = mean(magerr)
cont_std = stddev(mag)
jstart = min(mjd)
jend = max(mjd)
Larr = transpose(fltarr(nlc) + 1)
get_delta_tarr, mjd, delta_tarr = delta_tarr
Err = diag_matrix(magerr^2)
ii=indgen(nlc) # replicate(1,nlc)
jj=indgen(nlc) ## replicate(1,nlc)
diagarr = (ii eq jj)
inddiag = where(ii eq jj)
zydata = create_struct('mjd', mjd, 'mag', mag, 'magerr', magerr, $
                      'flux', flux, 'fluxerr', fluxerr, 'fluxivar', fluxivar, $
                      'nlc', nlc, 'Larr', Larr, 'delta_tarr', delta_tarr, $
											'Err', Err, 'inddiag', inddiag, 'tarr', mjd, 'y', mag)
end
