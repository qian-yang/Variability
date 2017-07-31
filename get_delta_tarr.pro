pro get_delta_tarr, tarr, delta_tarr = delta_tarr
  n = n_elements(tarr)
  ti = rebin(tarr, n, n, /sample)
  tj = transpose(ti)
  delta_tarr = abs(ti - tj)
  delta_tarr = double(delta_tarr)
  return
end
