function [R, Lags, RR]=maxxcorr2(comp,L,M)
%L... maximum lag
%M... length of data used for estimation
%comp... components
if nargin<3
    M=size(comp,2);
end

dim=size(comp,1);

RR=zeros(dim);
Lags=zeros(dim);
R=abs(comp*comp')/M;

for k=1:L
    Rk=abs(comp(:,1:M-k)*comp(:,k+1:M)')/(M-k);
    RR=RR+Rk;
    ind=find(R-Rk<0);
    R(ind)=Rk(ind);
    Lags(ind)=-k;    
    Rk=Rk';
    ind=find(R-Rk<0);
    R(ind)=Rk(ind);
    Lags(ind)=k;
end    
ind=find(eye(dim)>0); % vynulujeme diagonalu.
R(ind)=1;
