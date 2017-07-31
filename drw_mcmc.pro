;+
; Name: drw_mcmc
; Description: get the best fit parameters using MCMC (Metropolis-Hastings) method
; Operation are done by matrix, instead of circulation
; C = U^T U, where U is a upper triangular arrays calculated from LA_CHOLDC
; Verison: Qian Yang, July 26, 2017
;

pro drw_mcmc;, zydata, best_sigma = best_sigma, best_tau = best_tau
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
  file = 'r2.txt'
  read_data, file, zydata = zydata
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
  nwarker = 100
  log_sigma = 0.0
  log_tau = 0.0
  cov_sigma = 0.1
  cov_tau = 0.1
  result = replicate({log_sigma:0.d, log_tau:0.d, lnL:0.d, sigma:0.d, tau:0.d, prob:0.d, has:0}, nwarker)
  acceptance = 0.0
  kall = 0.0
  k = 0l
  stopped = 0
  while k le (nwarker - 1) do begin
    if (result[k].has eq 0) then begin
      result[k].has = 1
      result[k].log_sigma = log_sigma
      result[k].log_tau = log_tau
      result[k].sigma = 10.d^(result[k].log_sigma)
      result[k].tau = 10.d^(result[k].log_tau)
      sigma = result[k].sigma
      tau = result[k].tau
      ; DRW: sigma^2 * exp(-|ti - tj|/tau)
      C = (sigma)^2 * exp(-delta_tarr/(tau))
      C = C + Err
      lnlike_from_CL, y, C, L, inddiag, logL = logL
      result[k].lnL = logL
      result[k].prob = exp(result[k].lnL)
    endif
    ; mcmc
    prop_sigma = log_sigma + RANDOMN(seed) * cov_sigma
    prop_tau = log_tau + RANDOMN(seed) * cov_tau
    sigma = 10.d^prop_sigma
    tau = 10.d^prop_tau
    C = (sigma)^2 * exp(-delta_tarr/(tau))
    C = C + Err
    lnlike_from_CL, y, C, L, inddiag, logL = logL
    kall += 1.0
    if (logL ge result[k].lnL) then begin
      k += 1
      acceptance += 1
      log_sigma = prop_sigma
      log_tau = prop_tau
      print, log_sigma, log_tau
      stopped = 0
    endif else begin
      stopped += 1
      print, 'fail ', k, kall, acceptance/kall
      if (stopped gt 100) then break
    endelse
  endwhile
  print, 'kall ', kall
  print, 'acceptance ', acceptance
  print, 'acceptance ratio ', acceptance/kall
  tp = max(result.lnL, idbest)
  best_sigma_log = result[idbest].log_sigma
  best_tau_log = result[idbest].log_tau
  best_sigma = result[idbest].sigma
  best_tau = result[idbest].tau
  ;
  print, 'the best sigma and tau: ', best_sigma, best_tau

  filew = 'test6_mag_mcmc2.fits'
  spawn, 'rm -rf ' + filew
  mwrfits, result, filew
  print, SYSTIME()
end
