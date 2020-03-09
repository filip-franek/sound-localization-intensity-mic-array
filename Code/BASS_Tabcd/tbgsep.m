function [W W0 ISR]=tbgsep(x,L,thr)
%
% Blind separation for piecewise stationary sources
%
% Coded by Petr Tichavsky, October 2007. Updated June, 2009.
%
% Please cite
%
%  P. Tichavsky, A. Yeredor, "Fast Approximate Joint Diagonalization 
%      Incorporating Weight Matrices", IEEE Tr. on Signal Processing, 
%      vol. 57, no.3, pp. 878-891, March 2009. 
%
% x ..........  input signal
% L ..........  number of blocks
% W0 .........  initial separation by UWEDGE
% W ..........  final separation by WEDGE
% thr ........  threshold on variance of signals in blocks, 
%               below which the signal is considered to be too
%               low. Positive number between 0 and 1, close to zero.
%
[d N]=size(x);
if nargin<3
   thr=0.01;
end   
num_of_iterations = 3;
Nb=floor(N/L);
M0=zeros(d,d); M=[];
vars=zeros(1,L);
for l=1:L
    ll=(l-1)*Nb;
    Ma=x(:,ll+1:ll+Nb)*x(:,ll+1:ll+Nb)'/Nb;
    vars(1,l)=trace(Ma);
    M0=M0+Ma;
    M=[M Ma];
end
M0=M0/L;   
mez=thr*mean(vars)^2;
strong=find(vars.^2>mez);  %%% finds blocks with strong enough signals
lstrong=length(strong);
if lstrong<L & lstrong>1
   indstrong=reshape(ones(d,1)*(strong-1)*d+(1:d)'*ones(1,lstrong),1,d*lstrong);
   M=M(:,indstrong);
end   
[W0 Ms Rs]=uwedge2(M0,M);  
W=W0;
for in = 1:num_of_iterations
    [W Ms Rs]=wedge2(M0,M,Ms,Rs,W,5);
end
ISR=crb(Rs,Nb);
end %%%%%%%%%% of tbgsep

function [W Ms Rs]=uwedge2(M0,M,W0)
%
% an approximate joint diagonalization with uniform weights
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        M0 ... represents diagonalization constraint: 
%                W_est * M0 * W_est' has to have diagonal elements equal to 1.
%        W0 ... initial estimate of the demixing matrix, if available
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
improve=10;
if nargin<3
   [H E]=eig(M0);
   W=diag(1./sqrt(abs(diag(E))))*H';
else
   W=W0;
end  
Ms=M;  
Rs=zeros(d,L);
for k=1:L
      ini=(k-1)*d;
      M(:,ini+1:ini+d)=0.5*(M(:,ini+1:ini+d)+M(:,ini+1:ini+d)');
      Ms(:,ini+1:ini+d)=W*M(:,ini+1:ini+d)*W';
      Rs(:,k)=diag(Ms(:,ini+1:ini+d));
end 
crit=sum(Ms(:).^2)-sum(Rs(:).^2);  
C1=zeros(d,d); 
while improve>eps && iter<20
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
  critic=sum(Ms(:).^2)-sum(Rs(:).^2);
  improve=abs(critic-crit(end));
  crit=[crit critic]; 
  iter=iter+1;
end 
%crit
end %%%%%%%%%%%%%   of UWEDGE2 

function [W Ms Rs]=wedge2(M0,M,Ms,Rs,W,maxnumit)
%
% my approximate joint diagonalization with non-uniform weights
%
% Input: M .... the matrices to be diagonalized, stored as [M1 M2 ... ML]
%        M0 ... represents diagonalization constraint: 
%                W_est * M0 * W_est' has to have diagonal elements equal to 1.
%        W0 ... initial estimate of the demixing matrix, if available
%        maxnumit ... maximum number of iterations
% 
% Output: W .... estimated demixing matrix
%                    such that W_est * M_k * W_est' are roughly diagonal
%         Ms .... diagonalized matrices composed of W_est*M_k*W_est'
%         crit ... stores values of the diagonalization criterion at each 
%                  iteration
%
%
[d Md]=size(M);
L=floor(Md/d);
dd2=d*(d-1)/2;
Md=L*d;
iter=0;
eps=1e-7*dd2;
improve=1;
Rs0=Rs;
C1=zeros(d,d); 
while improve>eps && iter<maxnumit
  B=Rs.^2./Rs0*(1./Rs0)';
  aB=Rs./Rs0; B12=aB*aB';
  for id=1:d      
      C1(:,id)=Ms(:,id:d:Md)./Rs0*aB(id,:)';
  end
  D0=B.*B'-B12.^2;
  A0=eye(d)+(C1.*B-B12.*C1')./(D0+eye(d));
  W=A0\W;
  Raux=W*M0*W';
  aux=1./sqrt(diag(Raux));
  W=diag(aux)*W;  % normalize the result
  for k=1:L
     ini=(k-1)*d;
     Ms(:,ini+1:ini+d) = W*M(:,ini+1:ini+d)*W';
     Rs(:,k)=diag(Ms(:,ini+1:ini+d));
  end
  improve=sum(A0(:).^2)-d;
  iter=iter+1;
end  
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of wedge2

function ISR=crb(R,Nb)
%
[d L]=size(R);
ISR=zeros(d,d);
for k=1:d
    ISR(:,k)=sum(R./(ones(d,1)*R(k,:)),2);
end
ISR=ISR'./(ISR.*ISR'+eye(d)-L^2)/Nb;
ISR(eye(d)==1)=1;
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of crb
