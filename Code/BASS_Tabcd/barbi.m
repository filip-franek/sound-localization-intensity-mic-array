function [W W0 ISR]=barbi(x,L,ARo,rmax)
%
% Block-AutoRegressive Blind Identification (BARBI)
% =Blind separation for piecewise stationary AR sources
% =Block WASOBI separation
%
% Special cases: barbi(x,1,ARo) is another implementation of iwasobi(x,ARo)
%                barbi(x,L,0) is another implementation of BGL or bg_wedge(x,L)
%
% Coded by Petr Tichavsky, September 2008
%
% x ..........  input signal
% L ..........  number of blocks
% ARo ........  maximum AR order
% rmax .......  maximum modulus of AR poles 
% W0 .........  initial separation by UWEDGE
% W ..........  final separation 
%
[d N]=size(x);
Xmean=mean(x,2);
x=x-Xmean*ones(1,N); 
if nargin<4
   rmax=0.99;
end  
AR1=ARo+1;
num_of_iterations = 15;
eps=2e-6;
%Nb=floor((N-ARo)/L);   %%%% length of one block
fLN=floor(N/L);
Nb=floor(N/L)-ARo;
M0=zeros(d,d);       %%%% will be average covariance matrix with zero lag
M=zeros(d,d*L*(ARo+1));
ll=0; lm=0;
for l=1:L
    Ma=x(:,ll+1:ll+Nb)*x(:,ll+1:ll+Nb)'/Nb;
    M0=M0+Ma;
    M(:,lm+1:lm+d)=Ma;    
    for im=1:ARo
        lm=lm+d;
        Ma=x(:,ll+1:ll+Nb)*x(:,ll+im+1:ll+im+Nb)'/Nb;
        M(:,lm+1:lm+d)=0.5*(Ma+Ma');
    end    
    ll=ll+fLN;   
    lm=lm+d;
end
M0=M0/L;    
[W0 Rs Ms]=uwedge2(M0,M);  
W=W0;
Rx=zeros(AR1,d*L);
iter=0; crit=1;
while crit>eps && iter<num_of_iterations
  %  H = bwweights(R,L,rmax,eps0);
    for ir=1:AR1
        Rx(ir,:)=reshape(Rs(:,ir:AR1:end),d*L,1)';
    end  
    [ARC,sigmy]=armodel(Rx,rmax);  
    Rx3=ARC;
    for i=1:ARo
        Rx3(1:AR1-i,:)=Rx3(1:AR1-i,:)+(ones(AR1-i,1)*ARC(i+1,:)).*ARC(i+1:AR1,:);
    end
    Rx3(1,:)=Rx3(1,:)/2;
    Rx3=Rx3.*(ones(AR1,1)*(1./sigmy));
    Rx4=zeros(AR1*L,d);                 %% will be reshaped Rx3
    for il=1:L
        Rx4((il-1)*AR1+1:il*AR1,:)=Rx3(:,(il-1)*d+1:il*d);  
    end       
    B=Rs*Rx4;
    for id=1:d
      C1(:,id)=sum(Ms(:,id:d:end).*Rx4',2);
    end
    D0=B.*B'-diag(B)*diag(B)';
    A0=eye(d)+(C1.*B-diag(diag(B))*C1')./(D0+eye(d)); 
    W=A0\W;
    Raux=W*M0(:,1:d)*W';
    aux=1./sqrt(diag(Raux));
    W=diag(aux)*W;  % normalize the result
    for k=1:AR1*L
      ini=(k-1)*d;
      Ms(:,ini+1:ini+d) = W*M(:,ini+1:ini+d)*W';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
    end
    crit=(sum(abs(A0(:)))-d)/d^2;
    iter=iter+1;
  %  [iter crit]
end
ISR=CRLB6(ARC,sigmy,AR1,d)/Nb;
%t1 = [t1 cputime-time_start];
[isr poradi]=sort(sum(ISR,2)); % components are sorted according to their 
W=W(poradi,:);           % estimated isr (the most interesting ones go first)
W0=W0(poradi,:);
ISR=ISR(poradi,poradi);
%signals=W*x+(W*Xmean)*ones(1,N); %%% can be possible output as well
end %%%%%%%%%% of bwsep

function [W Rs Ms]=uwedge2(M0,M)
%
% an approximate joint diagonalization with uniform weights
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        M0 ... represents diagonalization constraint: 
%                W_est * M0 * W_est' has to have diagonal elements equal to 1.
% 
% Output: W  .... estimated demixing matrix
%                    such that W * M_k * W' are roughly diagonal
%         Ms .... diagonalized matrices composed of W*M_k*W'
%         crit ... stores values of the diagonalization criterion at each 
%                  iteration
% 
[d Md]=size(M);
L=floor(Md/d);
Md=L*d;
iter=0;
eps=1e-4;
[H E]=eig(M0);
W=diag(1./sqrt(abs(diag(E))))*H';
Ms=M;  
Rs=zeros(d,L);
for k=1:L
      ini=(k-1)*d;
      M(:,ini+1:ini+d)=0.5*(M(:,ini+1:ini+d)+M(:,ini+1:ini+d)');
      Ms(:,ini+1:ini+d)=W*M(:,ini+1:ini+d)*W';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
end 
crit=sum(Ms(:).^2)-sum(Rs(:).^2);  
while crit>eps && iter<100
  B=Rs*Rs';
  for id=1:d
      C1(:,id)=sum(Ms(:,id:d:Md).*Rs,2);
  end
  D0=B.*B'-diag(B)*diag(B)';
  A0=eye(d)+(C1.*B-diag(diag(B))*C1')./(D0+eye(d));
  W=A0\W;
  Raux=W*M0*W';
  aux=1./sqrt(abs(diag(Raux)));
  W=diag(aux)*W;  % normalize the result
  Ms=W*M;
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = Ms(:,ini+1:ini+d)*W';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
  crit=(sum(abs(A0(:)))-d)/d^2;
  iter=iter+1;
end 
end %%%%%%%%%%%%%   of UWEDGE2 
   
function [AR,sigmy]=armodel(R,rmax)
%
% to compute AR coefficients of the sources given covariance functions 
% but if the zeros have magnitude > rmax, the zeros are pushed back.
%
[M,d]=size(R);
switch M
    case 1
         sigmy=R; AR=ones(1,d);
    case 2
         v=R(2,:)./R(1,:);
         v=sign(v).*min([abs(v); rmax*ones(1,d)]);
         AR=[ones(1,d); -v];
         sigmy=R(1,:).*(1-v.^2);
    case 3
         den=1./(R(1,:).^2-R(2,:).^2);
         AR2=den.*R(2,:).*(R(3,:)-R(1,:)); 
         AR3=den.*(R(2,:).^2-R(1,:).*R(3,:));
         disc=sqrt(AR2.^2-4*AR3);
         v=[-AR2+disc; -AR2-disc]/2;
         vm=max(abs(v));
         ifac=rmax./max([vm; rmax*ones(1,d)]);
         AR=[ones(1,d); ifac.*AR2; ifac.^2.*AR3];   
         K=ones(1,d)-AR(3,:).^2;
         K1=-AR(2,:).*(ones(1,d)-AR(3,:))./K;
         sigmy=R(1,:).*K.*(1-K1.^2); 
    otherwise
AR=(levinson(R))'; %% can be replaced by line 171, if signal processing toolbox is not available
for id=1:d
%    AR(:,id)=[1; -toeplitz(R(1:M-1,id),R(1:M-1,id)')\R(2:M,id)];
    v=roots(AR(:,id)); %%% mimicks the matlab function "polystab"
    vs=0.5*(sign(abs(v)-1)+1);
    v=(1-vs).*v+vs./conj(v);
    vmax=max(abs(v));
    if vmax>rmax
       v=v*rmax/vmax;
    end   
    AR(:,id)=real(poly(v)'); %%% reconstructs back the covariance function
end 
Rs=ar2r(AR);
sigmy=R(1,:)./Rs(1,:);
end %%% of switch
end %%%%%%%%%%%%%%%%%%%%%%%  of armodel

function [ r ] = ar2r( a )
%%%%%
%%%%% Computes covariance function of AR processes from 
%%%%% the autoregressive coefficients using an inverse Schur algorithm 
%%%%% and an inverse Levinson algorithm (for one column it is equivalent to  
%%%%%      "rlevinson.m" in matlab)
% 
  if (size(a,1)==1)
      a=a'; % chci to jako sloupce
  end
  
  [p m] = size(a);    % pocet vektoru koef.AR modelu
  alfa = a;
  K=zeros(p,m);
  p = p-1;
  for n=p:-1:1
      K(n,:) = -a(n+1,:);
      for k=1:n-1
          alfa(k+1,:) = (a(k+1,:)+K(n,:).*a(n-k+1,:))./(1-K(n,:).^2);
      end
      a=alfa;
  end
%  
  r = zeros(p+1,m);
  r(1,:) = 1./prod(1-K.^2);
  f = r;
  b=f;
  for k=1:p 
      for n=k:-1:1
          K_n = K(n,:);
          f(n,:)=f(n+1,:)+K_n.*b(k-n+1,:);
          b(k-n+1,:)=-K_n.*f(n+1,:)+(1-K_n.^2).*b(k-n+1,:);
      end
      b(k+1,:)=f(1,:);
      r(k+1,:) = f(1,:);
  end       
end %%%%%%%%%%%%%%%%%%%%%%%%%%%  of ar2r

function ISR = CRLB6(ARC,sigmy,M,d);

%
% CRLB6(ARC) generates the CRLB for gain matrix elements (in term 
% of ISR) for blind separation of K Gaussian autoregressive sources 
% whose AR coefficients (of the length M, where M-1 is the AR order)
% are stored as columns in matrix ARC.

[M dL]=size(ARC);
L=floor(dL/d);

if M>1
   Rs=ar2r(ARC);
else
   Rs=ones(1,dL);
end   
power=Rs(1,:).*sigmy;
ISR=zeros(d,d); 
aveR=zeros(1,d);

for iL=1:L
    ill=(iL-1)*d+1:iL*d;
    aveR=aveR+power(1,ill);
    phi=zeros(d,d);
    for s=0:M-1
      for t=0:M-1
        phi=phi+(ARC(s+1,ill).*ARC(t+1,ill))'*Rs(abs(s-t)+1,ill);
      end
    end 
    ISR=ISR+diag(1./sigmy(1,ill))*phi*diag((sigmy(1,ill))); 
end
aveR=aveR/L;
denom=ISR'.*ISR+eye(d)-L^2; 
ISR=diag(aveR)*(ISR'./denom)*diag(1./aveR);;
ISR(eye(d)==1)=0;

end %%%%%%%%%%%%%%%%%%%%%%%%% of CRLB6


