;+
; NAME:
;   lnprob_func
; PURPOSE:
;   calculate the log-probabilities given parameters and data
;
; INPUT:
;   p: a vector of parameters
;   data: observational data, from read_data.pro
;
; OPTIONAL INPUTS:
;
; ARGS:
;   log: the parameters are in log space (alog10)
;   cont: if cont is set, there are only two parameters (sigma and tau)
;   scale: if scale is set, sigma = sigma_up * sqrt(tau/2.0)
;   lprior: if lprior is set, using lprior probability
;   pf: lprior format needed when lprior is set. ('Kelly09', 'Koz10', 'Zu11')
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;   logL: the log-probabilities
;
; COMMENTS:
;
; REVISION HISTORY:
;   29-Jul-2017  Written by Qian Yang, qianyang.astro@gmail.com
;
; Tudo:
;   add lprior probabiity
;-

function ln_lprior_cont, sigma, tau, pf, tt=tt, data=data, kargs=kargs
  my_neg_inf = -1.0e+300
  my_pos_inf = 1.0e+300
  tau_floor = 1.e-6
  tau_ceiling = 1.e+5
  sigma_floor = 1.e-6
  sigma_ceiling = 1.e+2
  logtau_floor = alog(tau_floor)
  logtau_ceiling = alog(tau_ceiling)
  logsigma_floor = alog(sigma_floor)
  logsigma_ceiling = alog(sigma_ceiling)
  lnpr = 0
  case pf of
  ; Kelly et al. 2009
  1 : lnpr -=  (tt)/tau
  ; Kozłlowski et al. 2010
  2 : lnpr -= alog(tau) + alog(sigma) ;sigma is sigma_up
  ; Zu et al. 2011 (JAVELIN)
  3 : begin
    lnpr -= alog(sigma); + alog(tau/tt)
    if (tau gt tt) then begin
        lnpr += - alog(tau/tt)
    endif else if (tau lt 0.001) then begin
      ; 86.4 seconds if input is in days
        lnpr += my_neg_inf
    endif else begin
        lnpr += - alog(tt/tau)
    endelse
    if check_args(kargs, 'taulimit') then begin
      taulimit = kargs.taulimit
      if (taulimit[1] le taulimit[0]) then taulimit = [data.tt, data.rj]
      if ((tau lt taulimit[0]) or (tau gt taulimit[1])) then begin
          lnpr += my_neg_inf
      endif
    endif
  end
  else : lnpr = -alog(tau)
  endcase
  return, lnpr
end

function ln_lprior_line, rj, lag, width, lagtobaseline=lagtobaseline, widtobaseline=widtobaseline
  my_neg_inf = -1.0e+300
  my_pos_inf = 1.0e+300
  if ~keyword_set(lagtobaseline) then lagtobaseline=0.3
  if ~keyword_set(widtobaseline) then widtobaseline=0.2
  lnpr = 0
  if (abs(lag) gt rj) then begin
    lnpr -= alog(abs(lag)/(lagtobaseline*rj))
    lnpr -= alog(abs(width)/(widtobaseline*rj))
  endif
  return, lnpr
end

function get_lnp, p, data, kargs=kargs
  ; unpack data
  tarr = data.tarr
  Larr = data.Larr
  y = data.y
  Err = data.Err
  inddiag = data.inddiag
  tt = data.tt
  rj = data.rj
  ; unpack params
  unpack_params, p, sigma=sigma, tau=tau, lag=lag, width=width, A=A, B=B, t1=t1, t2=t2, kargs=kargs
  get_S, tarr, Larr, sigma, tau, A=A, t1=t1, t2=t2, B=B, SS = SS
  C = SS + Err
  lnLike_from_CL, y, C, Larr, inddiag, logL = logL
  if check_args(kargs, 'lprior') then begin
    pf = kargs.pf
    lnpr = ln_lprior_cont(sigma, tau, pf, tt=tt, data=data, kargs=kargs)
    logL += lnpr
  endif
  if check_args(kargs, 'lp_line') then begin
    pf = kargs.pf
    lnpr = ln_lprior_line(rj, lag, width)
    logL += lnpr
  endif
  return, logL
end

pro lnprob_func, p, data, logL = logL, kargs=kargs
  ndim = n_elements(p[*, 0])
  n = n_elements(p[0, *])
  if (n eq 1) then begin
    logL = get_lnp(p, data, kargs=kargs)
  endif else begin
    logL = fltarr(n)
    for i = 0, n - 1 do begin
      logL[i] = get_lnp(p[*, i], data, kargs=kargs)
    endfor
  endelse
end
