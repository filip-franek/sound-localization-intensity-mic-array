function [order,clus,weights,unor,silh,Xdis] = RFCM2(R,c,m,init,e)

%Function for Relational Fuzzy C-means (output for convolutive demixing program, clusBGL)
%R - Dissimilarity matrix of data
%c - Number of clusters
%m - Fuzzyness weightening exponent
%e - Change of criterion

%checking the parameters given
%default parameters
if nargin<7, e = 1e-4; end
if nargin<4, init='wind'; end
if nargin<3, m = 1.2; end
if nargin<2, c = 2; end
if nargin<1, error('Data for clustering are missing'); end
ShowResults=1;
ShowHist=1;
OrderResults = 1;

dim = size(R,1);
maxit=100;
beta=0;
J=1;J0=0;                       % initialize stoping criterion
iter = 0;                       % iteration counter

% Similarity meassure into dissimilarity meassure
[X,Xbac]=sim2diss(R,'rc2');
%X=(R+R')/2;

%Show the histogram of the dissimilarity matrix
if(ShowHist)
   figure;subplot(1,3,1);hist(Xbac(:),20);title('Similarity matrix histogram:');
   subplot(1,3,2);hist(X(:),20);title('Dissimilarity matrix histogram:');
end   

% Initialize fuzzy partition matrix
f0=InitClus(c,dim,init);

% Iterate
while  (abs(J0-J) > e && iter<maxit)
   iter = iter + 1;
   J=J0;
   % Calculate centroids and square distance matrix (dsq)
   for lp = 1 : c,
      v(lp,:)=(f0(lp,:).^m)./sum(f0(lp,:).^m);
      dsq(lp,:)=(X*v(lp,:)'-((1/2)*v(lp,:)*X*v(lp,:)')*ones(dim,1))';
   end;
   %Spreading transformation (used in case of negative square distances)
   if(~isempty(find(dsq<0,1)))
      [dbeta,dsq]=CalcDBeta(v,dsq);
      beta=beta+dbeta;
      dbmat=ones(dim);dbmat=dbmat-diag(diag(dbmat));dbmat=dbmat*dbeta;
      X=X+dbmat;
   end
   %Update of fuzzy partition matrix
   dsqbac=dsq;
   dsq = (dsq+1e-10).^(-1/(m-1));
   f0 = (dsq ./ (ones(c,1)*sum(dsq)));

   %Calculate the criterion
   J0=sum(sum(f0.*dsqbac));

end
[vals,ind] = max(f0);
unor=ind;
[cluster,order] = sort(ind);
weights=f0(:,order);

%Formating data for clusBGL and check for empty clusters
unq=unique(ind);
ind=1;
for lp=unq
   tmp=find(cluster==lp);
   clus(ind,1) = tmp(1);
   clus(ind,2) = tmp(end)-tmp(1)+1;
   ind=ind+1;
end

%Calculate validation
silh=mean(Silhouette(unor,R));
%disp('Silhouette index:');disp(silh);

% save clustering variable
if OrderResults == 1
    X=R;
   X=(X+X')/2;
   Xdis=X(order,order);
end
%Image of Clustering
if(ShowResults)
   X=R;
   X=(X+X')/2;
   X=X(order,order);
   if(~ShowHist)
      figure;
      imagesc(X);
   else
      subplot(1,3,3);
      imagesc(X);      
   end
end
%%
%Function for calculation of spreading transformation value beta and
%respective update of dsq
function [dbeta,dsq2]=CalcDBeta(v,dsq)

[c,N]=size(v);
sqv=zeros(c,N);

for lp=1:c
   for lp2=1:N
      ek=zeros(1,N);ek(lp2)=1;
      tmp=v(lp,:)-ek;
      sqv(lp,lp2)=tmp*tmp';
   end
end

tmp=-2*(dsq)./(sqv);
dbeta=max(max(tmp));
dsq2=dsq+(dbeta/2)*sqv;

%%
%Function for different initialization types of relational clustering

function [f0]=InitClus(c,dim,type)

part=floor(dim/c);
f0=(0.5/c)*ones(c,dim);

switch (type)
   case 'mine'
      for lp=1:c
         for lp2=1:c
            if(lp2~=c)
               f0(lp,(lp2-1)*part+1:lp2*part)=mod((lp2-1)+(lp-1),c)/(c-1);
            else
               f0(lp,(lp2-1)*part+1:end)=mod((lp2-1)+(lp-1),c)/(c-1);
            end
         end
      end
   case 'mine2'
      tot=sum(1:c);mem=[1:c]*(1/tot);
      for lp=1:c
         for lp2=1:c
            if(lp2~=c)
               f0(lp,(lp2-1)*part+1:lp2*part)=mem(1+mod(lp+lp2-1,c));
            else
               f0(lp,(lp2-1)*part+1:end)=mem(1+mod(lp+lp2-1,c));
            end
         end
      end        
   case 'wind'
      for lp=1:c
         if(lp~=c)
            f0(lp,(lp-1)*part+1:lp*part)=0.5*(1+1/c);
         else
            f0(lp,(lp-1)*part+1:end)=0.5*(1+1/c);
         end
      end
end

%%
%Function for calculation of Silhouette Index
function silh=Silhouette(unor,R)

dim=size(R,1);
X=R;
%Symetrization
X=(X+X')/2;
%Similarity2dissimilarity
X=X/max(max(X));
X=X+diag(ones(1,dim));X=ones(dim)./X-diag(ones(1,dim));
X=X/max(max(X));

a=zeros(1,dim);b=zeros(1,dim);

for lp=1:dim
   same=find(unor==unor(lp));
   diff=find(unor~=unor(lp));
   a(lp)=sum(X(lp,same))/(length(same)-1);
   if(~isempty(diff))
      b(lp)=min(X(lp,diff));
   else
      b=0;
   end
end
silh=(b-a)./max(b,a);
%%
%Transformation of Similarity matrix into Dissimilarity matrix
function [X2,X]=sim2diss(R,kind)
   
   dim = size(R,1);
   grades=length(R(:)-dim)/2;

   X=R;
   %Data normalization
   X=(X+X')/2;

   %Similarity meassure into dissimilarity meassure
   X=X-min(min(X));X=X/max(max(X));
   switch(kind)
      case 'rec'
         X2=X+diag(ones(1,dim));X2=ones(dim)./X2;
      case 'rc2'
         b=0.3;
         X2=X+diag(ones(1,dim));X2=ones(dim)./(b+X2)-ones(dim)./(b+1);
      case 'lin'   
         X2=ones(dim)-X;
      case 'uni'   
         X=X+diag(100*ones(1,dim));
         Xs=sort(X(:),'descend');Xs=Xs(dim+1:end);Xl=length(Xs);X2=zeros(dim);
         for lp=1:grades
            tres=Xs(floor(lp*Xl/grades));
            X2(X<tres)=(1/grades)*lp;
         end
         X=X-diag(diag(X));
   end   
   X2=X2-diag(diag(X2));X2=X2-min(min(X2));X2=X2/max(max(X2));