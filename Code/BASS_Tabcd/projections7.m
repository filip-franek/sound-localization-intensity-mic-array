function R=projections7(comp,L,NFFT)
%L... number of GCC coefficients

if nargin<3
    NFFT=1024;
end

[dim N]=size(comp);

M=ceil(N/NFFT);

comp=[comp'; zeros(M*NFFT-N,dim)];

R=zeros(dim);

for k=1:M
    Fcomp=fft(comp((k-1)*NFFT+1:k*NFFT,:));
    Fcomp=sign(Fcomp);
    for i=1:dim
        for j=i+1:dim
            GCC=Fcomp(:,i).*conj(Fcomp(:,j));
            gcc=real(ifft(GCC));
            gcc=[gcc(1:L+1); gcc(end-L+1:end)];        
            %R(i,j,k)=sum(gcc.*sign(gcc));
            R(i,j)=R(i,j)+sum(gcc.*sign(gcc));
        end
    end
end
R=R+R';
