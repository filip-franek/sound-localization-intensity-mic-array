function [SBres perm permuted]=SBpermute3(SBres)
% Permutation of subbands of separated signals

verbose=true;

m = size(SBres,1); % number of microphones
N = size(SBres,2); % number of samples
n = size(SBres,3); % number of separated sources
nbands = size(SBres,4); % number of subbands

perm=zeros(nbands,n); % final permutations
permuted=zeros(1,nbands); % indicates if band was permuted

R=zeros(n*nbands);
for i=1:m
    R=R+abs(cov(abs(reshape(squeeze(SBres(i,:,:,:)),N,n*nbands))));
end
R=R.*~kron(eye(nbands),ones(n));

% Searching first two most matching bands and defining the reference
% permutation band
band1=1; % the first band is the reference one
perm(1,:)=1:n;
if verbose
    fprintf('Subband %d: %s\n',band1,num2str(1:n));
end
for i=1:nbands-1 
    val=max(triu(R((band1-1)*n+1:band1*n,:)));
    [val l]=max(val);
    band2=ceil(l/n); % the band to-be permuted
    p=greedyperm(R((band1-1)*n+1:band1*n,(band2-1)*n+1:band2*n));
    perm(band2,:)=p;
    if sum(p == 1:n)~=n
        permuted(band2)=1;
    end
    if verbose
        fprintf('Subband %d: %s\n',band2,num2str(p));
    end

    SBres(:,:,:,band2)=SBres(:,:,p,band2); %permute band2 and R correspondingly
    pom=R((band2-1)*n+1:band2*n,:);R((band2-1)*n+1:band2*n,:)=pom(p,:); 
    pom=R(:,(band2-1)*n+1:band2*n);R(:,(band2-1)*n+1:band2*n)=pom(:,p);
    % band1 won't be compared anymore
    R((band1-1)*n+1:band1*n,:)=0; 
    R(:,(band1-1)*n+1:band1*n)=0; 
    band1=band2;
end

%% Helping functions

function p=greedyperm(R)
n=size(R,1);
p=zeros(1,n);
for i=1:n
    [val indx]=max(R);
    [val l]=max(val);
    k=indx(l); % the maximum element is R(k,l)
    p(l)=k;
    R(k,:)=0;
    R(:,l)=0;
end