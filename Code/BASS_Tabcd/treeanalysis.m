function U=treeanalysis(g,x,M)
%Tree-structured Filter Bank analysis 

U=x;
for i=1:M
    Uold=U;
    [u0 u1]=qmfanalysis(g,Uold(1,:));
    U=zeros(2^i,length(u0));
    U(1,:)=u0;
    U(2,:)=u1;
    for j=2:2^(i-1)
        [u0 u1]=qmfanalysis(g,Uold(j,:));
        U(2*(j-1)+1,:)=u0;
        U(2*j,:)=u1;
    end
end
