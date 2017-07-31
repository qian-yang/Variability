pro call_func, p, logL = logL, data=data, kargs=kargs
  ;call_func, p, kargs=kargs, logL=logL, lg=lg, cont=cont, scale=scale, lprior=lprior, pf=pf
  ; lnprob_func, p, data, logL = logL
  ; lnprob_func, p, args, logL = logL
  lnprob_func, p, data, logL = logL, kargs=kargs;, lg=lg, cont=cont, scale=scale, lprior=lprior, pf=pf
end
