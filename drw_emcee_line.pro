pro drw_emcee_line
  print, SYSTIME()
  ; ======== box zoom input ========
  nwalkers = 100
  nchain = 100
  nburn = 100
  ; ===============================
  file1 = 'r.txt'
  file2 = 'g.txt'
  read_data_2, file1, file2, zydata = data, /ff ;flux
  ; mcmc
  ; sigma, tau, lag, width, A, B
  ndim = 6.0
  p0 = randomn(seed, ndim, nwalkers)
  kargs = create_struct('lg', 1, 'ln', 0, 'cont', 0, 'line', 1, 'scale', 0, 'lprior', 0, 'pf', 3, 'lp_line', 1, 'taulimit', [0, 0])
  emcee, p0, data=data, iterations=nchain, kargs=kargs, p=p, lnprob=lnprob
  emcee, p, data=data, iterations=nburn, kargs=kargs, lp0=lnprob, p=p, lnprob=lnprob, /save_chains, chains=chains, allp=allp
  result = replicate({sigma:0.d, tau:0.d, lag:0.d, width:0.d, A:0.d, B:0.d, lnL:0.d, log_sigma:0.d, log_tau:0.d}, nwalkers * nchain)
  if (kargs.lg) then begin
    result.log_sigma = transpose(chains[0, *])
    result.log_tau = transpose(chains[1, *])
    result.sigma = 10.d^result.log_sigma
    result.tau = 10.d^result.log_tau
  endif
  result.lag = transpose(chains[2, *])
  result.width = transpose(chains[3, *])
  result.A = transpose(chains[4, *])
  result.B = transpose(chains[5, *])
  result.lnL = allp

  filew = 'test_oneline.fits'
  spawn, 'rm -rf ' + filew
  mwrfits, result, filew
  print, SYSTIME()
end
