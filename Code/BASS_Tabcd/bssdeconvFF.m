function [signals, microphones]=bssdeconvFF(x,L,parametry)

%wtype.. 1..normal
%        2..constrained
%        3..binary
%        4..fuzzy
%        5..fuzzy constrained
%        6..fuzzy binary (the most likely cluster)

%clustype
%        'hclus'
%        'rfcm'

LengthOfICA=parametry.lengthofICA;
offset=parametry.offsetICA;
method=parametry.method;
wlevel=parametry.weighting;
wtype=parametry.wtype;
clustype=parametry.clustering;
mu=parametry.mu;
useiweights=parametry.useiweights;

m=size(x,1);

% fprintf('Length of data used for BSS: %d, offset: %d\n',LengthOfICA,offset);

% matrix defining the working subspace
dim=L*m;
N=length(x);
Lstar=(1+0.4*abs(mu-1)*log10(L))*L/mu;

X=zeros(dim,N); 
for j=1:m
    X((j-1)*L+1:j*L,:)=slaguerre(L,mu,x(j,:));
end;
y=X(:,offset+ceil(Lstar)-1:offset-1+LengthOfICA);

% waitbar(0,hwait,'ICA decomposition ...');
%ICA decomposition
if strcmp(method,'bgl')
   [We ISR]=tbgsep(y,floor(length(y)/300));
elseif strcmp(method,'efica')
   [We ISR]=efica(y,eye(size(y,1)));
elseif strcmp(method,'extefica')
   [We ISR]=befica(y,floor(length(y)/300),eye(size(y,1)));
elseif strcmp(method,'bwasobi')
   [We ISR]=barbi(y,floor(length(y)/300),1);
end
A=inv(We);
comp=We*y;

% waitbar(0.3,hwait,'Computing similarity of components ...');

%the distance matrix
if parametry.similarity==1
    R=projections4(comp,L);
elseif parametry.similarity==2
    R=projections7(comp,L);
elseif parametry.similarity==3
    R=projections8(comp,128);
end

%clustering
% waitbar(0.8,hwait,'Clustering ...');
if(strcmp(clustype,'hclus'))
   [p clst]=hclus(R,size(x,1));
   fprintf('Found %d clusters\n',size(clst,1))
elseif (strcmp(clustype,'rfcm'))
   [p clst weights_fuzzy]=RFCM2(R,size(x,1));
end


R=R(p,p);
ISR=ISR(p,p)-diag(diag(ISR));
We=We(p,:);
A=A(:,p);
signals=zeros(size(clst,1),size(x,2));

microphones=zeros(size(clst,1),length(x),size(x,1));

% waitbar(0.9,hwait,'Reconstruction ...');

if ~strcmp(method,'bgl')&&~strcmp(method,'bwasobi')
   iweights=sum(ISR,2);
   iweights=(1./iweights).^(1/2);
   iweights=iweights/max(iweights);
end
for i=1:size(clst,1) %forming a source from each cluster
   %forming all processed (delayed) sources
   incluster=clst(i,1):clst(i,1)+clst(i,2)-1;
   outcluster=setdiff(1:length(p),incluster);

   %Weights
   switch(wtype)
      case {1}     % Normal
         weights=(sum(R(:,incluster),2)./sum(R(:,outcluster),2));
         weights=(weights-min(weights)).^wlevel;
         if ~strcmp(method,'bgl')&&~strcmp(method,'bwasobi')&&useiweights
            weights=iweights.*weights;
         end
      case {2}
         if (strcmp(clustype,'hclus'))
             weights=compweights4d(comp(p,:),A,L,incluster,wlevel);
         else % RFCM
             weights=compweights4d(comp(p,:),A,L,find(weights_fuzzy(i,:)>0.4),wlevel);
         end
      case {3}       % Binary
         weights(incluster)=1;
         weights(outcluster)=0;
      case {4,5}     % Fuzzy
         %weights=(weights_fuzzy(i,:)./(1-weights_fuzzy(i,:))).^wlevel;         
         weights=weights_fuzzy(i,:).^wlevel;         
      case {6}       % Fuzzy Binary
         [tmp,ind]=max(weights_fuzzy);
         weights=zeros(1,size(weights_fuzzy,2));
         weights(ind==i)=1;         
   end
   weights=weights/max(weights);
   if (wtype==5)
      weights=AdjustWeights(weights);
   end

   H=A*diag(weights)*We;

   Xi=H*X;
   % microphones with separated source
   for j=1:m 
       microphones(j,:,i)=ilaguerre2(Xi((j-1)*L+1:j*L,:),mu);
   end

   %forming sources on microphones using the delay-and-sum beamformer
   [val j]=max(mean(microphones(:,:,i).^2,2));
   [val lag]=maxxcorr2(microphones(:,:,i),L);
   formmics=microphones(j,:,i);
   for k=[1:j-1 j+1:m]
       if lag(j,k)>0
          formmics=formmics+[zeros(1,lag(j,k)) microphones(k,1:end-lag(j,k),i)];
       else
          formmics=formmics+[microphones(k,-lag(j,k)+1:end,i) zeros(1,-lag(j,k))];
       end
   end
   formmics=formmics/m;


   %output signals
   signals(i,:)=formmics;
end





%% Laguerre - subspace
function U=slaguerre(L,mu,x)

N=length(x);

U=zeros(L,N);

U(1,:)=x;
U(2,:)=filter([0 mu],[1 mu-1],x); %low-pass

for i=3:L
    U(i,:)=filter([mu-1 1],[1 mu-1],U(i-1,:)); %all-pass
end

%% "inversion" of slaguerre
function x=ilaguerre2(U,mu)

[L N]=size(U);

Um=fliplr(U); %changing direction of time

x=Um(L,:);
for i=L:-1:3
    x=(Um(i-1,:)+filter([mu-1 1],[1 mu-1],x)); %non-causal inversion of an all-pass filter
end


x=fliplr(x);

x=filter([1/mu (mu-1)/mu],1,x); % inversion of the very first low-pass filter

x=([x(2:end) 0]+U(1,:))/L; 
