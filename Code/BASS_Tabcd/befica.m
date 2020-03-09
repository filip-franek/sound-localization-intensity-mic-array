function [W, ISR, icasig]=befica(X, num_seg, ini, SaddleTest, Uniform, Identity)
%% Block EFICA: [W, ISR, icasig]=befica(X, num_seg, ini, SaddleTest, Uniform, Identity)
% 
% Extended version of the EFICA algorithm designed for piecewise stationary non-Gaussian sources.
% version 1.02    release: 11.6.2009
%
% Input data:
% X ... mixed data dim x N, where dim is number of signals and N is number of
%       samples
% num_seg ... number of blocks in data
% ini ... starting point of the iterations
% SaddleTest ... if true (default), the test of saddle points on the demixed signals is done
% Uniform ... sets all weight parameters lambda to one: Uniform
%                Block EFICA (default is false)
% Identity ... if true (default), the identity function is used as the
%              third base function in the score function estimator
%
%Output data:
% W - demixing matrix produced by Block EFICA 
% ISR - ISR matrix estimator
% icasig - estimated independent components (normalized to unit scale)
%
% References
% [1] Z. Koldovský, J. Málek, P. Tichavský, Y. Deville, and S. Hosseini, 
%    "Blind Separation of Piecewise Stationary NonGaussian Sources", accepted for publication 
%     in Signal Processing, 2009.
%
% [2] Z. Koldovský, J. Málek, P. Tichavský, Y. Deville, and S. Hosseini, 
%    "Extension of EFICA Algorithm for Blind Separation of Piecewise Stationary Non Gaussian Sources",
%    ICASSP 2008, Las Vegas, April 2008.

[dim N]=size(X);

%% Default settings
if nargin<6
    Identity=true;
end
if nargin<5
    Uniform=false;
end
if nargin<4
    SaddleTest=true;
end
if nargin<3
    ini=eye(dim);
end
if nargin<2
    num_seg=1;          %Number of segments of signals 
end

g='rat1';               %Contrast function used for the symmetric part of Block EFICA (EEF1)
epsilon=0.0001;         %Stop criterion
fineepsilon=1e-7;       %Stop criterion for finetunings (EEF2)
MaxIt=100;              %Maximum number of iterations for symmetric part of Block EFICA (EEF1)
MaxItAfterSaddleTest=30;%Maximum number of iterations after a saddle point was indicated
FinetuneMaxIt=50;       %Maximum number of finetuning iterations (EEF2)
min_correlation=0.9;    %additive noise...0.75, noise free... 0.95, turn off (unstable)...0
test_of_saddle_points_nonln='rat2'; %nonlinearity used for the test of saddle points
global num_basis;  
if Identity
    num_basis=3;            %Number of basis functions considered in the score function estimator
else
    num_basis=2;
end
seg_len=fix(N/num_seg); %The length of segments 
PerfEst=true;           %If true, perform score function estimation in every iteration (EEF2)
%minseglen=100;          %Minimal allowed segment length 
extraseg=rem(N,num_seg);%The length of the last segment

% if (N/num_seg<minseglen)
%     error('Input data are too short for selected number of segments!'); 
% end    

fprintf('Starting Block EFICA, dim=%d, N=%d, M=%d, ',dim,N,num_seg);

W = symdecor(ini);
NumIt=0;
TotalIt=0;
crit=zeros(1,dim);
repeat=1;               %Further interations - symmetric part

%% Removing mean of data
Xmean=mean(X,2);
X=X-Xmean*ones(1,N);

%% Preprocessing of data
C = cov(X');
CC = C^(-1/2);
Z = CC*X;       %Preprocessed data

%% Symmetrix FastICA with the Test of saddle points 
G=zeros(dim,N); 
EGS=zeros(dim,num_seg);
while repeat
    while (1-min(crit)>epsilon && NumIt<MaxIt)
        NumIt=NumIt+1;  
        Wold=W;
        switch g
            case 'estm'
                U=W'*Z;
                for j=1:dim
                    [EGS(j,:) G(j,:)]=ScoreEstim(U(j,:),num_seg,true,zeros(num_basis,num_seg),Identity); 
                end
                W=Z*G'/N-ones(dim,1)*mean(EGS,2)'.*W;
            case 'tanh'
                hypTan = tanh(Z'*W);
                W=Z*hypTan/N-ones(dim,1)*sum(1-hypTan.^2).*W/N;
            case 'pow3'
                W=(Z*((Z'*W).^ 3))/N-3*W;
            case 'rat1'
                U=Z'*W;
                Usquared=U.^2;
                RR=4./(4+Usquared);
                Rati=U.*RR;
                Rati2=Rati.^2;
                dRati=RR-Rati2/2;  
                nu=mean(dRati);
                hlp=Z*Rati/N;
                W=hlp-ones(dim,1)*nu.*W;
            case 'rat2'
                U=Z'*W;
                Ua=1+sign(U).*U;
                r1=U./Ua;
                r2=r1.*sign(r1);
                Rati=r1.*(2-r2);
                dRati=(2./Ua).*(1-r2.*(2-r2));  
                nu=mean(dRati);
                hlp=Z*Rati/N;
                W=hlp-ones(dim,1)*nu.*W;
            case 'gaus'
                U=Z'*W;
                Usquared=U.^2;
                ex=exp(-Usquared/2);
                gauss=U.*ex;
                dGauss=(1-Usquared).*ex;
                W=Z*gauss/N-ones(dim,1)*sum(dGauss).*W/N;
        end
        TotalIt=TotalIt+dim;        
        W=symdecor(W);
        crit=abs(sum(W.*Wold));
    end %while iteration
    
    if repeat==1
        fprintf('Iterations: %d\n',NumIt);
    elseif repeat==2
        fprintf('   Test of saddle points positive: %d iterations\n',NumIt);
    end
    
    repeat=0; %Do not repeat the Test of saddle points anymore
    
    %%%The Test of saddle points of the separated components
    if SaddleTest
        SaddleTest=false; %%The SaddleTest may be done only one times
        u=Z'*W;
        switch test_of_saddle_points_nonln
            case 'tanh'
                table1=(mean(log(cosh(u)))-0.37456).^2;
            case 'gaus'
                table1=(mean(ex)-1/sqrt(2)).^2;
            case 'rat1'
                table1=(mean(2*log(4+u.^2))-3.1601).^2;
            case 'rat2'
                table1=(mean(u.^2./(1+sign(u).*u))-0.4125).^2;
            case 'pow3'
                table1=(mean((pwr(u,4)))-3).^2;
        end
         %  applying the round-Robin tournament scheme for parallel processing 
        dimhalf=floor((dim+1)/2); dim2=2*dimhalf; 
        da=[1:dim2 2:dim2];      %%% auxiliary array    
        for delay = 0:dim2-2
          ii=[1 da(dim2-delay+1:3*dimhalf-delay-1)];
          jj=da(4*dimhalf-delay-1:-1:3*dimhalf-delay);
          if dim2>dim
             i0=dimhalf-abs(delay-dimhalf+1/2)+1/2;
             ii(i0)=[]; jj(i0)=[];  % the pair containing index dim2 must be deleted
          end 
         ctrl0=table1(ii)+table1(jj);
         z1=(u(:,ii)+u(:,jj))/sqrt(2);
         z2=(u(:,ii)-u(:,jj))/sqrt(2);
         switch test_of_saddle_points_nonln
                case 'tanh'
                     ctrl=(mean(log(cosh(z1)))-0.37456).^2 ...
                          +(mean(log(cosh(z2)))-0.37456).^2;
                case 'gaus'
                     ctrl=(mean(exp(-z1.^2/2)-1/sqrt(2))).^2 ...
                          +(mean(exp(-z2.^2/2)-1/sqrt(2))).^2;
                case 'rat1'
                     ctrl=(mean(2*log(4+z1.^2))-3.1601).^2 ...  
                          +(mean(2*log(4+z2.^2))-3.1601).^2;  
                case 'rat2'
                     ctrl=(mean(z1.^2./(1+sign(z1).*z1))-0.4125).^2 ...
                          +(mean(z2.^2./(1+sign(z2).*z2))-0.4125).^2;
                case 'pow3'
                     ctrl=(mean((pwr(z1,4)))-3).^2 ...
                          +(mean((pwr(z2,4)))-3).^2;
         end
         indexes=find(ctrl>ctrl0);
         if length(indexes)>0
            irot=ii(indexes); jrot=jj(indexes);
              %bad extrems indicated
           % fprintf('  Block EFICA: rotating components: %d\n', [irot jrot]);
            u(:,irot)=z1(:,indexes); 
            u(:,jrot)=z2(:,indexes); 
            Waux=W(:,irot);
            W(:,irot)=(W(:,irot)+W(:,jrot))/sqrt(2);
            W(:,jrot)=(Waux-W(:,jrot))/sqrt(2);
            NumIt=0;MaxIt=MaxItAfterSaddleTest;%status=1;
            repeat=2; %continue in iterating - the test of saddle points is positive
         end %if length(indeces)>0
        end% for delay
    end %if SaddleTest
    crit=zeros(1,dim);
end %while repeat


%Estimated signals by the Symmetric FastICA with the test of saddle points (EEF1)
Wsymm=W;

s=W'*Z;   

[mu nu beta]=params(reshape(s(:,1:(num_seg-1)*seg_len)',seg_len,(num_seg-1)*dim)',g);
[tmpmu tmpnu tmpbeta]=params(s(:,(num_seg-1)*seg_len+1:N),g);           %Parameters of last segment with possible extraseg samples
mu=reshape(mu,num_seg-1,dim)';nu=reshape(nu,num_seg-1,dim)';beta=reshape(beta,num_seg-1,dim)';     
mu(:,num_seg)=tmpmu;nu(:,num_seg)=tmpnu;beta(:,num_seg)=tmpbeta;     

%% One-unit FastICA finetuning for piecewise stationary signals 
% utilizes Pham's score function estimator

lambda=ones(dim,num_seg);  
for j=1:dim
    w=W(:,j);
    estpars=zeros(num_basis,num_seg);%%% initialization of the score function estimator
    estimate=true;
    wold=zeros(dim,1);
    nit=0;
    while ((1-abs(w'*wold)>fineepsilon) && (nit<FinetuneMaxIt) &&...
            (abs(W(:,j)'*w)>min_correlation))
        [EGS G estpars]=ScoreEstim((w'*Z),num_seg,estimate,estpars,Identity);
        tmpEXG=(w'*Z).*G;
        EXG(1:num_seg-1)=mean(reshape(tmpEXG(1:(num_seg-1)*seg_len),seg_len,num_seg-1));
        EXG(num_seg)=mean(tmpEXG((num_seg-1)*seg_len+1:N));                 %The mean of last whole segment with added samples of the extraseg(length(extraseg)<seg_len)
        if(~PerfEst)
            estimate=false; 
        end    
        if(~Uniform)
            lambda(j,:)=(EGS-EXG)*inv(num_seg*diag(EGS)-EXG'*EXG);
            lambdas=ones(seg_len,1)*lambda(j,:);
            lambdas=lambdas(:)';
            if(extraseg~=0)                                                 %Lambda vector in case that there are samples in extraseg
                lambdas=[lambdas ones(1,extraseg)*lambda(j,num_seg)];    
            end    
            wold=w;          
            w=Z*(G.*lambdas)'/N-mean(EGS.*lambda(j,:))*w;            
        else
            wold=w;          
            w=Z*G'/N-mean(EGS)*w;            
        end    
        w=w/norm(w);
        nit=nit+1;
        TotalIt=TotalIt+1;
    end
    if abs(W(:,j)'*w)>min_correlation,
        %%%%  the signal has not been changed too much, so the finetuning was likely successful                
        [beta(j,:) G]=ScoreEstim(w'*Z,num_seg,false,estpars,Identity);
        nu(j,:)=beta(j,:);
        tmpmu=(w'*Z).*G;
        mu(j,1:num_seg-1)=mean(reshape(tmpmu(1:(num_seg-1)*seg_len),seg_len,num_seg-1));
        mu(j,num_seg)=mean(tmpmu((num_seg-1)*seg_len+1:N));   %The mean of last whole segment with added samples of the extraseg(length(extraseg)<seg_len)
        W(:,j)=w;
    else
        lambda(j,:)=ones(1,num_seg);   
    end
end

fprintf('Total number of iterations/component: %.1f\n',TotalIt/dim);

%% ISR of one-unit de-mixing vector estimates 
ISR1U=((mean(lambda.^2.*beta,2)-mean(lambda.*mu,2).^2)./(mean(lambda.*(nu-mu),2).^2))*ones(1,dim);

%% Refinement (EEF3)
Werr=zeros(dim); 
for k=1:dim
    ccc=ISR1U(k,:)./(ISR1U(:,k)'+1);ccc(k)=1;
    WW=W*diag(ccc);
    sloupce=setdiff(1:dim,find(sum(WW.^2)/max(sum(WW.^2))<1e-7)); % removing almost zero rows
    M=WW(:,sloupce);
    M=symdecor(M);
    if sum(sloupce==k)==1
        Wextefica(:,k)=M(:,sloupce==k);
    else
        w=null(M');
        if size(w,2)==1
            Wextefica(:,k)=w(:,1);
        else % there are probably more than two gaussian components => use the old result to get a regular matrix
            Wextefica(:,k)=Wsymm(:,k);
        end
    end
    %Estimate variance of elements of the gain matrix
    Werr(k,:)=(ISR1U(k,:).*(ISR1U(:,k)'+1))./(ISR1U(k,:)+ISR1U(:,k)'+1);
end
Werr=Werr-diag(diag(Werr));

%% Output
ISR=Werr/N;
W=Wextefica'*CC;
%Wsymm=Wsymm'*CC;
icasig=W*X+(W*Xmean)*ones(1,N);



%% Fast power computation
function x=pwr(a,n)
x=a;
for i=2:n
    x=x.*a;
end

%% Fast symmetric orthogonalization
function W=symdecor(M)
[V D]=eig(M'*M);
W=M*(V.*(ones(size(M,2),1)*(1./sqrt(diag(D)'))))*V';

%% Estimation of the moments mu, nu, and beta
function [mu,nu,beta]=params(s,g)
[dim N]=size(s);
global num_basis;
switch g
    case 'estm'
        G=zeros(dim,N);
        EGS=zeros(dim,1);
        Identity=evalin('caller','Identity');
        for j=1:dim
            [EGS(j,:) G(j,:)]=ScoreEstim(s(j,:),1,1,zeros(num_basis,1),Identity); 
        end
        mu=mean(s.*G,2);
        nu=EGS;
        beta=EGS;
    case 'tanh'
        mu=mean(s.*tanh(s),2);
        nu=mean(1./cosh(s).^2,2);
        beta=mean(tanh(s).^2,2);
    case 'rat1'
        ssquare=s.^2;
        mu=mean(ssquare./(1+ssquare/4),2);
        nu=mean((1-ssquare/4)./(ssquare/4+1).^2,2);
        beta=mean(ssquare./(1+ssquare/4).^2,2);
    case 'rat2'
        r1=s./(1+s.*sign(s));
        r2=r1.*sign(r1);
        Rati=r1.*(2-r2);
        dRati=(2./(1+s.*sign(s))).*(1-r2.*(2-r2));  
        mu=mean(s.*Rati,2);
        nu=mean(dRati,2);
        beta=mean(Rati.^2,2);
    case 'gaus'
        aexp=exp(-s.^2/2);
        mu=mean(s.^2.*aexp,2); 
        nu=mean((1-s.^2).*aexp,2); 
        beta=mean((s.*aexp).^2,2);
    case 'pow3'
        mu=mean(s.^4,2); 
        nu=3*ones(dim,1); 
        beta=mean(s.^6,2);
end

%% Score function estimator
function [EGS,psi,estpars] = ScoreEstim(data,num_seg,doest,estpars,Identity)
datlen=length(data);
seg_len=fix(datlen/num_seg);
EGS=zeros(1,num_seg);
psi=zeros(1,datlen);

%computation of the nonlinearities x^3 and x/(1+6*|x|)^2 and their
%derivatives
data=data';
x1pow2=data.*data; 
x1pow3=x1pow2.*data;
x1pow4=x1pow2.*x1pow2;
r1=1./(1+6*data.*sign(data)); r2=r1.*r1;
x2exp=data.*r2;
x2expp=r2-12*x2exp.*sign(x2exp).*r1;
if Identity
    h=zeros(datlen,7);
    h(:,1)=x1pow2;
    h(:,2)=x1pow4;
    h(:,3)=x1pow2.*r2;
    h(:,4)=x1pow3.^2;
    h(:,5)=x1pow4.*r2;
    h(:,6)=x2exp.^2;
    h(:,7)=x2expp;
    HH=reshape(sum(reshape(h(1:seg_len*(num_seg-1),:),seg_len,(num_seg-1)*7)),num_seg-1,7); 
    HH(end+1,:)=sum(h(seg_len*(num_seg-1)+1:end,:));
else
    h=[x1pow3.^2 x1pow4.*r2 x2exp.^2 x1pow2 x2expp];
    HH=reshape(sum(reshape(h(1:seg_len*(num_seg-1),:),seg_len,(num_seg-1)*5)),num_seg-1,5); 
    HH(end+1,:)=sum(h(seg_len*(num_seg-1)+1:end,:));
end

for i=1:num_seg  
    first=(i-1)*seg_len+1;
    if(i==num_seg)                  %Last segment could have different length than the others
        last=datlen;
        seg_len=datlen-first+1;
    else    
        last=i*seg_len;
    end
    H=HH(i,:);
    if Identity
        dt=H(1)*(H(4)*H(6)-H(5)^2)-H(2)*(H(2)*H(6)-2*H(3)*H(5))-H(3)^2*H(4);
        a11 = H(6)*H(4) - H(5)^2;
        a12 = H(3)*H(5)-H(2)*H(6);
        a13 = H(2)*H(5)-H(3)*H(4);
        a22 = H(1)*H(6)-H(3)^2;
        a23 = H(2)*H(3)-H(1)*H(5);
        a33 = H(1)*H(4) - H(2)^2;
        w(1)=(seg_len*a11+3*H(1)*a12+a13*H(7))/dt;
        w(2)=(seg_len*a12+3*H(1)*a22+a23*H(7))/dt;
        w(3)=(seg_len*a13+3*H(1)*a23+a33*H(7))/dt;
        EGS(i)=(seg_len*w(1)+3*H(1)*w(2)+H(7)*w(3))/seg_len; %Note that, here, E[psi']=E[psi^2] !
    else
        hhT=[H([1 2]);H([2 3])];
        sumdh=[3*H(4); H(5)]; 
        invhhT=inv(hhT);
        w=invhhT*sumdh;
        EGS(i)=w'*sumdh/seg_len; %Note that, here, E[psi']=E[psi^2] !
    end    
    if doest  %if true, re-estimate optimum linear combination of functions
        estpars(:,i)=w;
    else                
        w=estpars(:,i);        
    end     

    if Identity
        psi(first:last)=w(1)*data(first:last)+w(2)*x1pow3(first:last)+w(3)*x2exp(first:last);
    else
        psi(first:last)=w(1)*x1pow3(first:last)+w(2)*x2exp(first:last);
    end
end



