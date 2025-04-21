%simulation，10 times
clear
[N,A,kgd1,kgd2,edges,tri] = dataReal2();

rng(10)
Omega=randn(1,N);
mean(kgd1)
mean(kgd2)

tmax=2000;
dt=0.01;
eqstart1=1;

tic

k1=[-60:-1];
k2=[-60:-1];

cnt=1;
for i=1:length(k1)
    for j=1:length(k2)
        kk(cnt,:)=[k1(i),k2(j)];
        cnt=cnt+1;
    end
end
cnt=cnt-1
kks=kk;

kk1=kks(:,1);
kk2=kks(:,2);
Ck1 = parallel.pool.Constant(kk1);
Ck2 = parallel.pool.Constant(kk2);

% tic
% [R1,RR1,theta] = kurad2h(tmax,kgd1,kgd2,edges,tri,Omega,...
%     -30,-10,dt,eqstart1);   
% toc
parpool(80)
parfor i=1:cnt
    [R1(i),RR1(i),theta(i,:)] = kurad2(tmax,kgd1,kgd2,edges,tri,Omega,...
        Ck1.Value(i),Ck2.Value(i),dt,eqstart1);    
end

toc

