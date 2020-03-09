function [g0, g1, h0, h1]=qmf(g)
% design of quadrature mirror filter bank
N=length(g);
g0=g;
g1=(-1).^(0:N-1)' .* g;
h0=2*g0;
h1=-2*g1;