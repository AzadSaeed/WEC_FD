function [a_, b_, fe_] = MS_calculations(Omega,pdata)

n_wec  = pdata.General.n_wec;
x      = pdata.WEC.x;
y      = pdata.WEC.y;
nw     = length(Omega);
Radius = pdata.WEC.Radius;
Draft  = pdata.WEC.Draft;
depth  = pdata.General.depth;
g      = pdata.General.g;
rho    = pdata.General.rho;
w      = pdata.WEC.w;

[a_, b_, fe_]  = MyMSCoefficients(n_wec,[x,y],nw,Radius,Draft,depth,g,rho,w,pdata.MS);


end