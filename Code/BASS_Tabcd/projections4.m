function R=projections4(comp,L)
%L... number of delays

[dim N]=size(comp);

R=zeros(dim);
comp2=[zeros(dim,L) comp zeros(dim,L)];

cQ=[];
for k=-L:L
    cQ=[cQ comp*comp2(:,L+k+1:N+L+k)'];
end

Fcomp=fft([comp'; zeros(2^nextpow2(N)-N,dim)]);
%r=real(ifft(Fcomp.*conj(Fcomp)));
r=real(length(Fcomp)*zoomfft(Fcomp.*conj(Fcomp),length(Fcomp),2^nextpow2(2*L+1),0,1));
QQTinvALL=FTI(r(1:2*L+1,:),2*L,dim);


for k=1:dim
%    c=comp(k,L+1:end);
%    for l=1:2*L+1, Q(l,l:l+N-1)=c;end
%    Q=datafordeconvolution(comp(k,L+1:end),2*L);
%    Q=Q(:,1:size(comp,2));
%    QQT=toeplitz(r(1:2*L+1,k)'); %%% zbytecne
    cmps=[1:k-1 k+1:dim];
%    A=(comp(cmps,:)*Q(:,1:size(comp,2))')*inv(QQT);
%    A=(comp(cmps,:)*Q(:,1:size(comp,2))')*QQTinvALL(:,(2*L+1)*(k-1)+1:(2*L+1)*k);
%    Au1=comp(cmps,:)*Q(:,1:size(comp,2))';
%    raux=real(ifft((conj(Fcomp(:,k))*ones(1,dim)).*Fcomp));
    Au2=fliplr(cQ(cmps,k:dim:end));
%    Au2=[raux(end-L+1:end,cmps)' raux(1:L+1,cmps)'];
    A=Au2*QQTinvALL(:,(2*L+1)*(k-1)+1:(2*L+1)*k);
%    E=mean((A*Q).^2,2);
%    E2=diag(A*QQT*A')/(N+2*L);
    %%% ekvivalentne E2=diag(A*Au2')/(N+2*L); nebo E2=sum(A.*Au2,2)/(N+2*L);
%    real(ifft(fft(A',length(Fcomp)).*(Fcomp(:,k)*ones(1,dim-1))));
%    R(k,cmps)=E2';
     R(k,cmps)=sum(A.*Au2,2)'/(N+2*L);
end    

function Q3L=FTI(r,M,d)
%
% implements Fast Toeplitz inversion as
%  Q1=[];
%  for k=1:d
%      Q1=[Q1 inv(toeplitz(r(:,k)))];
%  end
%    
Md=(M+1)*d;
[A,E]=levinson(r);    %%%%%%%%%%%% needs signal processing toolbox
Q3L=zeros(M+1,d*(M+1));
Q3L(1,1:M+1:Md)=1;
Q3=ones(d,1);
for m=1:M
        Q3(:,1:m+1)=[zeros(d,1) Q3(:,1:m)]+(A(:,m+1)*ones(1,m+1)).*A(:,1:m+1)...
                    -(A(:,M-m+2)*ones(1,m+1)).*[zeros(d,1) A(:,M+1:-1:M-m+2)];
        Q3a=[Q3 zeros(d,M-m)]';
        Q3L(m+1,:)=Q3a(:)';
end
for k=1:d
        Q3a=Q3L(:,(k-1)*(M+1)+1:k*(M+1));
        Q3L(:,(k-1)*(M+1)+1:k*(M+1))=(Q3a+Q3a'-diag(diag(Q3a)))/E(k);
end    


