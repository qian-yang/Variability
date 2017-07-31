;+
; Name: get_S
; Description: similar to spear_covfunc.f90 in JAVELIN
; Operation are done by matrix, instead of circulation
; C = U^T U, where U is a upper triangular arrays calculated from LA_CHOLDC
; Verison: Qian Yang, July 26, 2017
;

function check_symm, S, name
  ST = transpose(S)
  del = S - ST
  res = 0
  if ((min(del) eq 0.0) and (max(del) eq 0.0)) then begin
    print, name + ' is symmetric.'
    res = 1
  endif else begin
    print, name + ' is not symmetric!!!'
    res = -1
  endelse
  return, res
end

function get_delta_tarr, tarr1, tarr2, S = S
  n1 = n_elements(tarr1)
  n2 = n_elements(tarr2)
  S = dblarr(n1, n2)
  ti = rebin(tarr1, n1, n2, /sample)
  tj = transpose(rebin(tarr2, n2, n1, /sample))
  delta_tarr = ti - tj
  return, delta_tarr
end

function cont_cont, tarr1, tarr2, sigma, tau
  ; tarr1 cont
  ; tarr2 cont
  delta_t = get_delta_tarr(tarr1, tarr2)
  S = sigma^2 * exp(-abs(delta_t)/tau)
  return, S
end

function cont_line, tarr1, tarr2, A, t1, t2, sigma, tau
  ; tarr1 cont
  ; tarr2 cont+line
  if (t1>t2) then begin
    tp = t1
    t1 = t2
    t2 = tp
  endif
  delta_t = get_delta_tarr(tarr1, tarr2, S = S)
  tL = delta_t - t2
  tH = delta_t - t1
  exp_tL_tau = exp(tL/tau)
  exp_tH_tau = exp(tH/tau)
  ind1 = where(tL gt 0, nind1)
  if (nind1 gt 0) then S[ind1] = 1.0/exp_tL_tau[ind1] - 1.0/exp_tH_tau[ind1]
  ind2 = where(tH lt 0, nind2)
  if (nind2 gt 0) then S[ind2] = exp_tH_tau[ind2] - exp_tL_tau[ind2]
  ind3 = where((tL le 0) and (tH ge 0), nind3)
  if (nind3 gt 0) then S[ind3] = 2.0 - exp_tL_tau[ind3] - 1.0/exp_tH_tau[ind3]
  S = S * tau * sigma^2 * A / (t2 - t1)
  return, S
end

; function line_cont, S
;   S2 = transpose(S)
;   return, S2
; end

function line1_line2, tarr1, tarr2, A, t1, t2, B, t3, t4, sigma, tau
  ; tarr1 cont+line1
  ; tarr2 cont+line2
  ; t1>t2
  if (t1 gt t2) then begin
    tp = t1
    t1 = t2
    t2 = tp
  endif
  ; t3>t4
  if (t3 gt t4) then begin
    tp = t3
    t3 = t4
    t4 = tp
  endif
  ; (t2 - t1) > (t4 - t3)
  if ((t4 - t3) gt (t2 - t1)) then begin
    tp1 = t1
    tp2 = t2
    tp = A
    t1 = t3
    t2 = t4
    A = B
    t3 = tp1
    t4 = tp2
    B = tp
  endif
  delta_t = get_delta_tarr(tarr1, tarr2, S = S)
  tL  = delta_t - (t2-t3)
  tM1 = delta_t - (t2-t4)
  tM2 = delta_t - (t1-t3)
  tH  = delta_t - (t1-t4)
  ; tH > tM2 > tM1 > tL
  exp_tL_tau = exp(-abs(tL)/tau)
  exp_tH_tau = exp(-abs(tH)/tau)
  exp_tM1_tau = exp(-abs(tM1)/tau)
  exp_tM2_tau = exp(-abs(tM2)/tau)
  S = exp_tL_tau + exp_tH_tau - exp_tM1_tau - exp_tM2_tau
  ind1 = where((tH gt 0) and (tM2 le 0), nind1)
  if (nind1 gt 0) then S[ind1] += 2.0 * tH[ind1]/tau
  ind2 = where((tM2 gt 0) and (tM1 le 0), nind2)
  if (nind2 gt 0) then S[ind2] += 2.0 * (t4 - t3)/tau
  ind3 = where((tM1 gt 0) and (tL lt 0), nind3)
  if (nind3 gt 0) then S[ind3] += -2.0 * tL[ind3]/tau
  S = S * tau^2 * sigma^2 * A * B / (t2 - t1) / (t4 - t3)
  return, S
end

function line_line, tarr1, tarr2, A, t1, t2, sigma, tau
  ; tarr1 cont+line1
  ; tarr2 cont+line2
  ; t1>t2
  if (t1 gt t2) then begin
    tp = t1
    t1 = t2
    t2 = tp
  endif
  ; t3 = t1, t4 = t2, B = A
  delta_t = get_delta_tarr(tarr1, tarr2, S = S)
  tL  = delta_t - (t2-t1)
  tM = delta_t
  tH  = delta_t - (t1-t2)
  ; tH > tM > tL
  exp_tL_tau = exp(-abs(tL)/tau)
  exp_tH_tau = exp(-abs(tH)/tau)
  exp_tM_tau = exp(-abs(tM)/tau)
  S = exp_tL_tau + exp_tH_tau - 2.0 * exp_tM_tau
  ind1 = where((tH gt 0) and (tM lt 0), nind1)
  if (nind1 gt 0) then S[ind1] += 2.0 * tH[ind1]/tau
  ind3 = where((tM ge 0) and (tL lt 0), nind3)
  if (nind3 gt 0) then S[ind3] += -2.0 * tL[ind3]/tau
  S = S * tau^2 * sigma^2 * A^2 / (t2 - t1)^2
  return, S
end

pro get_S, tarr, Larr, A, t1, t2, sigma, tau, B, SS = SS
  ind1 = where(Larr[0, *] eq 1)
  tarr1 = tarr[ind1]
  ind2 = where(Larr[1, *] eq 1)
  tarr2 = tarr[ind2]
  S11 = cont_cont(tarr1, tarr1, sigma, tau)
  S12 = B * cont_cont(tarr1, tarr2, sigma, tau) + cont_line(tarr1, tarr2, A, t1, t2, sigma, tau)
  S21 = transpose(S12)
  S22_cl = cont_line(tarr2, tarr2, A, t1, t2, sigma, tau)
  S22_lc = transpose(S22_cl)
  S22_ll = line_line(tarr2, tarr2, A, t1, t2, sigma, tau)
  S22 = B^2 * cont_cont(tarr2, tarr2, sigma, tau) + B * (S22_cl + S22_lc) + S22_ll
  S1 = [S11, S12]
  S2 = [S21, S22]
  SS = transpose([transpose(S1), transpose(S2)])
  return
end
