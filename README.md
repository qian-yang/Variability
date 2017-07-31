Variability
====
This is an IDL package to analysis QSO variability and do photometric reverberation mapping (photo-RM) <br>

Notes
----
The core code are learnt from the JAVELIN code (please see https://bitbucket.org/nye17/javelin) and the Python emcee package (http://dan.iel.fm/emcee/current/). There is a modified version of JAVELIN code, including some features such as, adding structure function and changing the transfer functions by Haowen (please see https://github.com/pkuzhw387/Modified_JAVELIN)<br>

Current Version
----
Qian Yang <br>
July 30, 2017 <br>

Already done (July 30, 2017)
----
Damped random walk (DRW) model QSO continuum variability (probability calculated in lnlike_from_CL.pro) <br>
Continuum + emission line with lag model (covariance matrix in get_S.pro) <br>
Maximum the probability, using zoom in boxes (run drw_zoomin.pro) <br>
MCMC using Metropolis-Hastings method (run drw_mcmc.pro) <br>
MCMC using method learn from emcee (run drw_emcee.pro) <br>
MCMC (emcee.pro debug) <br>

TODO (July 30, 2017) <br>
----
Structure function (compare with previous works) <br>
Continuum lag (no transfer function) <br>
DRW in mag, while scale in flux <br>
Consider host galaxy contribution (not change), which is different from emission line contribution (has lag and transfer function) <br>
Predict light curves (simulate LSST light curve) <br>
Try different format of transfer functions <br>
DRW using different bands<br>
Two emission lines lag<br>
Multi-band analysis <br>
Define conditions (S/N, time baseline, time interval, line ratio) <br>

Contract
-----
Any question please feel free to email Qian (qianyang.astro@gmail.com).
