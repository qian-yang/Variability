function mag2flux, mag
	; flux = (10.0^((-0.4) * mag)) * (3.631*(1e-20)) ; erg s-1 Hz-1 cm-2
	; flux = (10.0^((-0.4) * mag)) * (3.631e6) ; mJy (Jy Jansky)
	flux = (10.0^((-0.4) * mag)) * (3.631e3) ; Jy (Jy Jansky)
  return, flux
end

function magerr2fluxerr, mag, magerr, flux
	fluxerr = flux * 0.4 * alog(10.d) * magerr
	return, fluxerr
end

pro read_data_2, file1, file2, zydata = zydata, ff=ff
readcol, file1, format = 'd, d, d', mjd1, mag1, magerr1
readcol, file2, format = 'd, d, d', mjd2, mag2, magerr2
;
ind1 = where((mjd1 gt 0) and (mag1 gt 0) and (magerr1 gt 0), nlc1)
mjd1 = mjd1[ind1]
mag1 = mag1[ind1]
magerr1 = magerr1[ind1]
sort_index = sort(mjd1)
tp = mjd1[sort_index]
tt = median(abs(tp[1:nlc1-1] - tp[0:nlc1-2]))
jstart = min(mjd1)
jend = max(mjd1)
rj = jend - jstart
;
ind2 = where((mjd2 gt 0) and (mag2 gt 0) and (magerr2 gt 0), nlc2)
mjd2 = mjd2[ind2]
mag2 = mag2[ind2]
magerr2 = magerr2[ind2]
;
flux1 = mag2flux(mag1)
fluxerr1 = magerr2fluxerr(mag1, magerr1, flux1)
fluxivar1 = 1.0/(fluxerr1^2)
;
flux2 = mag2flux(mag2)
fluxerr2 = magerr2fluxerr(mag2, magerr2, flux2)
fluxivar2 = 1.0/(fluxerr2^2)
;
mjd = [mjd1, mjd2]
mag = transpose([mag1, mag2])
magerr = [magerr1, magerr2]
flux = transpose([flux1, flux2])
fluxerr = [fluxerr1, fluxerr2]
print, median(flux), median(fluxerr)
Larr = [transpose([fltarr(nlc1) + 1, fltarr(nlc2) + 0]), transpose([fltarr(nlc1) + 0, fltarr(nlc2) + 1])]
nlc = nlc1 + nlc2
;
yy = mag
Err = diag_matrix(magerr^2)
if keyword_set(ff) then begin
	yy = flux
	Err = diag_matrix(fluxerr^2)
endif
ii=indgen(nlc) # replicate(1,nlc)
jj=indgen(nlc) ## replicate(1,nlc)
diagarr = (ii eq jj)
inddiag = where(ii eq jj)
writearr, [transpose(mjd1), transpose(flux1), transpose(fluxerr1)], 'RMID589_r_flux.txt'
writearr, [transpose(mjd2), transpose(flux2), transpose(fluxerr2)], 'RMID589_g_flux.txt'
;
zydata = create_struct( $
				'nband', 2, $
				'mjd1', mjd1, 'mag1', mag1, 'magerr1', magerr1, $
				'flux1', flux1, 'fluxerr1', fluxerr1, 'nlc1', nlc1, $
				'mjd2', mjd2, 'mag2', mag2, 'magerr2', magerr2, $
				'flux2', flux2, 'fluxerr2', fluxerr2, 'nlc2', nlc2, $
				'mjd', mjd, 'mag', mag, 'magerr', magerr, $
				'flux', flux, 'fluxerr', fluxerr, $
				'nlc', nlc, 'Larr', Larr, $
				'tarr', mjd, 'y', yy, 'Err', Err, 'inddiag', inddiag, $
				'tt', tt, 'rj', rj)
end
