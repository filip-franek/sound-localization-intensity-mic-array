function [u0 u1]=qmfanalysis(g,x)

[g0 g1]=qmf(g);
u0=filter(g0,1,x);
u1=filter(g1,1,x);
u0=u0(1:2:length(u0));
u1=u1(1:2:length(u1));

