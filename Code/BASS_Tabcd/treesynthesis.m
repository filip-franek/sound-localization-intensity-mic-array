function U=treesynthesis(g,U)
%Tree-structured Filter Bank synthesis

M=log2(size(U,1));
for i=M:-1:1
    Uold=U;
    [v0 v1]=qmfsynthesis(g,Uold(1,:),Uold(2,:));
    U=zeros(2^(i-1),length(v0));
    U(1,:)=v0+v1;
    for j=2:2^(i-1)
        [v0 v1]=qmfsynthesis(g,Uold(2*(j-1)+1,:),Uold(2*j,:));
        U(j,:)=v0+v1;
    end
end
