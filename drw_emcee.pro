pro drw_emcee
  print, SYSTIME()
  ; ======== input ========
  nwalkers = 50
  nchain = 100
  nburn = 100
  ; ===============================
  file = 'r2.txt'
  read_data, file, zydata = data
  ; mcmc
  ndim = 2.0
  p0 = randomn(seed, ndim, nwalkers)
  kargs = create_struct('lg', 1, 'ln', 0, 'cont', 1, 'line', 0, 'scale', 0, 'lprior', 1, 'pf', 3, 'taulimit', [0, 0])
  emcee, p0, data=data, iterations=nchain, kargs=kargs, p=p, lnprob=lnprob
  emcee, p, data=data, iterations=nburn, kargs=kargs, lp0=lnprob, p=p, lnprob=lnprob, /save_chains, chains=chains, allp=allp
  result = replicate({sigma:0.d, tau:0.d, lnL:0.d, log_sigma:0.d, log_tau:0.d, ln_sigma:0.d, ln_tau:0.d}, nwalkers * nchain)
  if kargs.ln then begin
    result.ln_sigma = transpose(chains[0, *])
    result.ln_tau = transpose(chains[1, *])
    result.log_sigma = result.ln_sigma/alog(10.d)
    result.log_tau = result.ln_tau/alog(10.d)
  endif
  if (kargs.lg) then begin
    result.log_sigma = transpose(chains[0, *])
    result.log_tau = transpose(chains[1, *])
    result.ln_sigma = result.log_sigma * alog(10.d)
    result.ln_tau = result.log_tau * alog(10.d)
  endif
  result.sigma = 10.d^result.log_sigma
  result.tau = 10.d^result.log_tau
  result.lnL = allp

filew = 'test_lg_prior.fits'
spawn, 'rm -rf ' + filew
mwrfits, result, filew
print, SYSTIME()
end
