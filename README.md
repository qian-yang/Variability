Variability
====
This is an IDL package to analysis QSO variability and do photometric reverberation mapping (photo-RM) <br>

Version
----
Qian Yang <br>
July 25, 2017 <br>

References
----
The core code are learnt from the JAVELIN code (please see https://bitbucket.org/nye17/javelin) and the Python emcee package (http://dan.iel.fm/emcee/current/). There is a modified version of JAVELIN code, including some features such as, adding structure function and changing the transfer functions by Haowen (please see https://github.com/pkuzhw387/Modified_JAVELIN)<br>

Already done (July 30, 2017)
----
Damped random walk (DRW) model QSO continuum variability (core in lnlike_from_CL.pro) <br>
Continuum + emission line with lag model (core in get_S.pro) <br>
Maximum the probability, using zoom in boxes (in drw_zoomin.pro) <br>
MCMC using Metropolis-Hastings method (in drw_mcmc.pro) <br>
MCMC using method learn from emcee (in emcee.pro) <br>
MCMC (emcee.pro debug) <br>

TODO (July 30, 2017) <br>
----
Structure function (compare with previous works) <br>
Statistic the lag between different continuum bands <br>
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