function w=compweights4d(C,A,L,incluster,penalty)
% C ... independent components
% A ... mixing matrix
% L ... the length of demixing filter (length of blocks)
% incluster ... indices of components that should be taken into account

M=size(C,1);
m = M/L;

CC=zeros(M,M,2*L-1);
for ell=0:L-1
    CC(:,:,ell+L)=C(:,ell+1:end-3*L+ell)*C(:,1:end-3*L)';
    CC(:,:,L-ell)=CC(:,:,ell+L)';
end

H=zeros(M);
for ell=1:M
    for ell2=1:M
        for r=1:m    
            for k=1:L
                for p=1:L
                    H(ell,ell2)=H(ell,ell2)+A((r-1)*L+p,ell)*A((r-1)*L+k,ell2)*CC(ell,ell2,p-k+L);
                end
            end
        end
    end
end

H=(A'*A).*CC(:,:,L)-H/L;

%H=diag(1./diag(H))*H;

w=zeros(M,1);

HH=H(incluster,incluster);
M2=size(HH,1);
HH=HH/norm(HH)+(1/penalty)*(eye(M2)-(1/M2)*ones(M2));
[V D]=eig(HH);
[~,ind]=min(diag(D));
% Prvky V(:,ind) by mely mit vsechny stejne znamenko, aby to mohly
% byt vahy. Je zajimave a dobre, ze tomu tak skoro vzdycky je. Pritom
% rigorozni zduvodneni zatim nezname.
w(incluster)=abs(V(:,ind)); 
w=w/max(w); % normovani
