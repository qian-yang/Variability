;+
; NAME:
;   unpack_params
;
; PURPOSE:
;   unpack parameters
;
; INPUT:
;   p: a vector of parameters
;
; OPTIONAL INPUTS:
;
; KEYWORDS:
;   log: the parameters are in log10 space (alog10)
;   cont: if cont is set, there are only two parameters (sigma and tau)
;   scale: if scale is set, sigma = sigma_up * sqrt(tau/2.0)
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;   sigma: a long-term standard deviation of variability
;   tau: a damping timescale
;   A: the flux ratio of line to continuum if there is a band with cont+line
;   t1: t1 for top-hat transfer function
;   t2: t2 for top-hat transfer function
;   B: the flux ratio between line band (cont+line) and cont band
;
; COMMENTS:
;
; REVISION HISTORY:
;   29-Jul-2017  Written by Qian Yang, qianyang.astro@gmail.com
;
; Tudo:
;-

pro unpack_params, p, sigma=sigma, tau=tau, lag=lag, width=width, A=A, B=B, t1=t1, t2=t2, kargs=kargs
  ; log
  if check_args(kargs, 'lg') then begin
    sigma = 10.d^p[0]
    tau = 10.d^p[1]
    ; p = 10.d^p
  endif
  if check_args(kargs, 'ln') then begin
    sigma = exp(p[0])
    tau = exp(p[1])
    ; p = exp(p)
  endif
  if check_args(kargs, 'line') then begin
    lag = p[2]
    width = p[3]
    A = p[4]
    B = p[5]
    t1 = lag - width/2.d
    t2 = lag + width/2.d
  end
  if check_args(kargs, 'scale') then begin
    sigma = sigma * sqrt(tau/2.d)
  endif
end
