function R=projections8(comp,NFFT)
% Similarity of independent component based on their coherence

if nargin<2
    NFFT=512;
end

[dim N]=size(comp);

M=floor(N/NFFT);
%comp=[comp'; zeros(M*NFFT-N,dim)];

Fcomp=zeros(NFFT/2+1,M,dim);
Fn=zeros(NFFT/2+1,dim);
for i=1:dim
    Fcomp(:,:,i)=spectrogram(comp(i,:),hanning(NFFT),0,NFFT);
    %Fn(:,i)=sqrt(sum(Fcomp(:,:,i).*conj(Fcomp(:,:,i)),2));
    Fn(:,i)=sum(Fcomp(:,:,i).*conj(Fcomp(:,:,i)),2);
end

R=zeros(dim);

for i=1:dim
    for j=i+1:dim
        coherence=sum(Fcomp(:,:,i).*conj(Fcomp(:,:,j)),2);
        %coherence=sqrt(coherence.*conj(coherence))./Fn(:,i)./Fn(:,j);%./...
        coherence=(coherence.*conj(coherence))./Fn(:,i)./Fn(:,j);%./...
        R(i,j)=mean(coherence);
        %R(i,j)=mean(coherence.*hamming(length(coherence)));
    end
end
R=R+R';
