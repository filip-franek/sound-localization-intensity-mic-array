function [v0 v1]=qmfsynthesis(g,u0,u1)

[g0 g1 h0 h1]=qmf(g);

v0=reshape([u0;zeros(1,length(u0))],1,2*length(u0));
v1=reshape([u1;zeros(1,length(u1))],1,2*length(u1));
v0=filter(h0,1,v0);
v1=filter(h1,1,v1);

