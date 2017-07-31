pro drw_zoomin, zydata, best_sigma = best_sigma, best_tau = best_tau
  ; ======== box zoom input ========
  maxiter = 4
  ns = 21
  nt = 42
  step_sigma0 = 0.5
  step_tau0 = 0.5
  bin_sigma0 = 5
  bin_tau0 = 10
  ; ===============================
  print, SYSTIME()
  nlc = zydata.nlc
  tarr = zydata.mjd
  y = zydata.mag
  L = zydata.Larr
  magerr = zydata.magerr
  Err = diag_matrix(magerr^2)
  get_delta_tarr, tarr, delta_tarr = delta_tarr
  ii=indgen(nlc) # replicate(1,nlc)
  jj=indgen(nlc) ## replicate(1,nlc)
  diagarr = (ii eq jj)
  inddiag = where(ii eq jj)
  all = ns * nt * maxiter
  print, all
  result = replicate({log_sigma:0.d, log_tau:0.d, lnL:0.d, zoomin:-1, sigma:0.d, tau:0.d, prob:0.d}, all)
  best_sigma_log = 0.0
  best_tau_log = 0.0
  for iter = 0, maxiter - 1 do begin
    step_sigma = step_sigma0/(10.d^iter)
    step_tau = step_tau0/(10.d^iter)
    bin_sigma = bin_sigma0/(10.d^iter)
    bin_tau = bin_tau0/(10.d^iter)
    log_sigma_arr = best_sigma_log + dindgen(ns) * step_sigma - bin_sigma ; -11 11
    log_tau_arr = best_tau_log + dindgen(nt) * step_tau - bin_tau ; -11 11
    know = ns * nt * iter
    for i = 0, ns - 1 do begin
      for j = 0, nt - 1 do begin
        k = know + i * nt + j
        result[k].log_sigma = log_sigma_arr[i]
        result[k].log_tau = log_tau_arr[j]
        result[k].sigma = 10.d^(result[k].log_sigma)
        result[k].tau = 10.d^(result[k].log_tau)
        result[k].zoomin = iter
        sigma = result[k].sigma
        tau = result[k].tau
        ; DRW: sigma^2 * exp(-|ti - tj|/tau)
        C = (sigma)^2 * exp(-delta_tarr/(tau))
        C = C + Err
        lnlike_from_CL, y, C, L, inddiag, logL = logL
        result[k].lnL = logL
      endfor
    endfor
    index = where(result.zoomin eq iter)
    tp = max(result[index].lnL, idbest)
    best_sigma_log = result[index[idbest]].log_sigma
    best_tau_log = result[index[idbest]].log_tau
  endfor
  best_sigma = result[index[idbest]].sigma
  best_tau = result[index[idbest]].tau
  result.prob = exp(result.lnL)
  ;
  print, 'the best sigma and tau: ', best_sigma, best_tau

  filew = 'test6_mag_zoomin5_convolve.fits'
  spawn, 'rm -rf ' + filew
  mwrfits, result, filew
  print, SYSTIME()
end
