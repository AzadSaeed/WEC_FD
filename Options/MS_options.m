function pdata = MS_options(pdata)

pdata.MS.maxs  = 5;                          % order of interaction
pdata.MS.maxm  = 5;                          % order of series truncation for main fluid
pdata.MS.maxn  = 40;                         % order of series truncation for fluid below body
pdata.MS.maxj  = 15;                         % order of series truncation for diffraction problem
pdata.MS.fmaxw = pdata.General.wMax;         % max frequency in rad/s (for wave.m)
pdata.MS.fresw = pdata.General.nbfreq;       % number of frequencies for wave.m evaluations (linearly spaced from 0 to fmaxw)
pdata.MS.fres  = 50;                         % number of frequencies for spectral evaluations (linearly spaced from 0 to fmaxw)

end