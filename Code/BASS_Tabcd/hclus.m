function [order,clus_out,qof] = hclus(R,ncls)
% Agglomerative hierarchical clustering with averaged linkage strategy

d = size(R,1);
R=R+R';

%R=R+R';

% distance matrix
D = -.5*R;
%D = -.5*(max(R,R'));

% initial clusters are the individual components
cluster = zeros(d,d);
cluster(1,:) = 1:d;

% merge clusters until there is no more clusters to merge
cdim=ones(1,d);
for i = 2:d,
    % initialize the clusters
    cluster(i,:) = cluster(i-1,:);
%   slow computation of the distances    
%     for j=1:d
%         for k=1:d
%             DDD(j,k)=mean(mean(D(cluster(i,:)==cluster(i,j),cluster(i,:)==cluster(i,k))));
%         end
%     end
    % find pair of clusters that are closest to each other
    % prefer merging of smaller clusters
    [mval,mrow] = min(D./min(cdim'*ones(1,d),ones(d,1)*cdim));
    [mval,mcol] = min(mval);
    mrow = mrow(mcol);
    % merge the clusters
    C1 = find(cluster(i,:)==cluster(i,mrow));
    C2 = find(cluster(i,:)==cluster(i,mcol));
    cindex = min(cluster(i,[C1 C2]));
    cluster(i,[C1 C2]) = cindex;
    %dimension of the new cluster
    d1 = length(C1);d2 = length(C2);dn = d1+d2;    
    cdim([C1 C2]) = dn; 
    % Make sure that we will not try to merge them again
    D(C1,C2) = 0;
    D(C2,C1) = 0;
    % compute the new distances
    D(:,[C1 C2]) = ((d1*D(:,C1(1))+d2*D(:,C2(1)))/dn)*ones(1,dn);
    D([C1 C2],:) = ones(dn,1)*((d1*D(C1(1),:)+d2*D(C2(1),:))/dn);
end

% select the best partition level - estimate the number of clusters

qof=zeros(1,d);
for i = d-ncls+1:d-1,
    cindex = unique(cluster(i,:));
    cq = [];
    cqo = [];
    for j = 1:length(cindex)
        Cin = find(cluster(i,:)==cindex(j));
        Cout = find(cluster(i,:)~=cindex(j)); %setdiff(1:d,Cin);
        %evaluation of intra-cluster dependencies
%         RR=R(Cin,Cin);
%         [val p]=max(triu(RR,1));
%         [val q]=max(val);        
%         cq(j)=0;
%         for k=1:size(RR,1)-1
%             RR(:,q)=0;
%             [val q]=max(RR(q,:));
% %             RR(:,q)=0;
%             cq(j)=cq(j)+val;
% % %            imagesc(RR)
%         end
%         cq(j)=cq(j)+max(RR(:,q));
%         cq(j)=cq(j)/size(RR,1);
        cq(j) = mean(mean(R(Cin,Cin))); %sum(max(R(Cin,Cin),[],2));
        
        %evaluation of inter-cluster dependencies
        cqo(j) = mean(mean(R(Cin,Cout),2));
%         end
    end
%    warning off MATLAB:divideByZero
    qof(i) = mean(cq./cqo);
%    warning on MATLAB:divideByZero
end
[qof,plevel] = max(qof);

%plevel=d-ncls+1;

% ordering of components and setting up the output format
D = -.5*(R);
D = D-eye(size(D));
[cluster,order] = sort(cluster(plevel,:));
cindex = unique(cluster);
D = D(order,order);
for i = 1:length(cindex)
    tmp = find(cluster==cindex(i));
    clus_out(i,1) = tmp(1);
    clus_out(i,2) = tmp(end)-tmp(1)+1;
    Cin = find(cluster==cindex(i));
    Cout = setdiff(1:d,Cin);
    if length(Cin) == 1,
        num = -1;
        den = min(min(D(Cin,Cout)));
        dlen = 1;
        clus_out(i,3) = sum(sum(R(order(Cin),order(Cout))))*(d-1)/(dlen*(d-dlen));
    elseif isempty(Cout),
        den = -1;
        num = max(max(D(Cin,Cin)));
        clus_out(i,3) = NaN;
    else
        num = max(max(D(Cin,Cin)));
        den = min(min(D(Cin,Cout)));
        dlen = length(Cin);
        clus_out(i,3) = sum(sum(R(order(Cin),order(Cout))))*(d-1)/(dlen*(d-dlen));
    end
    clus_out(i,4) = num/den;
end