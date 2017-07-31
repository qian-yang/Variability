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

function tri_determ, Matrix, inddiag
  ;similar to chodet_from_tri in cholesky_utils.py
  U = Matrix
  LA_CHOLDC, U, STATUS = test
  U_part = U[inddiag]
  det_C_ln = 2.0 * total(alog(U_part))
  return, det_C_ln
end

;+
; Name: lnlike_from_CL
; Description: similar to _lnlike_from_U in lcmodel.py in Javelin
; C = U^T U, where U is a upper triangular arrays calculated from LA_CHOLDC
; Verison: Qian Yang, July 25, 2017
;
pro lnlike_from_CL, y, C, L, inddiag, logL = logL
  ; |C|
  det_C = LA_DETERM(C, /CHECK)
  if ((det_C le 0.0) or FINITE(det_C, /INFINITY)) then begin
    det_C_ln = tri_determ(C, inddiag)
  endif else begin
    ; ln|C|
    det_C_ln = alog(det_C) ; ==> the first item
  endelse
  ; C don't have to be a diag matrix
  Ccopy = C
  LA_CHOLDC, Ccopy, STATUS = test
  if (test eq 0) then begin
    ; C a = y ==> a = C^-1 y
    a = LA_CHOLSOL(Ccopy, y)
    ; C b = L ==> b = C^-1 L
    b = LA_CHOLSOL(Ccopy, L)
    ; Cp = L^T C^-1 L = Cq^-1 = L^T b
    Lr = transpose(L)
    Cp = Lr ## b
    n = n_elements(Cp)
    if (n eq 1) then begin
      Cp = Cp[0]
      det_Cp = Cp
      d = (1.0/Cp) * Lr
      det_Cp_ln = alog(det_Cp) ; ==> the second item
    endif else if (n gt 1) then begin
      det_Cp = LA_DETERM(Cp, /CHECK)
      if ((det_Cp le 0.0) or FINITE(det_Cp, /INFINITY)) then begin
        det_Cp_ln = tri_determ(Cp, inddiag)
      endif else begin
        det_Cp_ln = alog(det_Cp)
      endelse
      ; Cp d = L^T ==> d = Cp^-1 L^T
      Cp_copy = Cp
      LA_CHOLDC, Cp_copy, STATUS = test
      d = LA_CHOLSOL(Cp_copy, Lr)
    endif else begin
      print, "Warning: n_elements(Cp) < 0"
    endelse
    ; e = C^-1 L Cq L^T C^-1 y = (C^-1 L) (Cq L^T) (C^-1 y) = b d a
    e = b ## d ## a
    ; f = C^-1 y - C^-1 L C_p^-1 L^T C^-1 y = C_v^-1 y
    f = a - e
    ; h = y^T C_v^-1 y
    h = transpose(y) ## f
    ; logL
    logL = -0.5d * det_C_ln - 0.5d * det_Cp_ln - 0.5d * h
    ; print, logL
  endif else begin
    logL = -999.0
  endelse
  if (FINITE(logL, /INFINITY) or FINITE(logL, /NAN)) then logL = -999.0
  return
end
