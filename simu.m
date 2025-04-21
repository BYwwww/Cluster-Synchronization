clear;
N=1000;
tmax=2000;     
m=20;
NX=100;
a=[-9:0];
x=logspace(3,1,21);
k2=[-x,a,1:30];
k1=[-30:30];
dt=0.01;
eqstart1=1;


tic
parpool(80)

load('net.mat')


%-------------------------------------
disp(['begin..'])
len1=length(k1);
len2=length(k2);
cnt=1;
for i=1:len1
    for j=1:len2
        kk(cnt,:)=[k1(i),k2(j)];
        cnt=cnt+1;
    end
end
cnt=cnt-1



kk1=kk(:,1);
kk2=kk(:,2);
Ck1 = parallel.pool.Constant(kk1);
Ck2 = parallel.pool.Constant(kk2);

%alpha=1;
%[kgd210,tri10] = graphpp(N,A,kgd2,tri,alpha);

% M=sum(kgd1)/2
% ntris=sum(kgd2)/3
% [Ar,kgr,kg2r,edgesr,trinr] = graphpnm(A,N,M,ntris,alpha)

Ckg1 = parallel.pool.Constant(kgd1);
Ckg2 = parallel.pool.Constant(kgd2);
Cedges = parallel.pool.Constant(edges);
Ctri = parallel.pool.Constant(tri);


% Ckg210 = parallel.pool.Constant(kgd210);
% Ctri10 = parallel.pool.Constant(tri10);


parfor i=1:cnt

    [R1(i),RR1(i),theta(i,:)] = kurad2(tmax,kgd1,kgd2,edges,tri,...
        Omega,Ck1.Value(i),Ck2.Value(i),dt,eqstart1);   
    %[R110(i),RR110(i)] = kurad2(tmax,Ckg1.Value,Ckg210.Value,Cedges.Value,Ctri10.Value,...
       % Omega,Ck1.Value(i),Ck2.Value(i),dt,eqstart1);  
end

toc


