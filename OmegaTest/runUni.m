clear;
N=1000;
tmax=1000;     %时长
m=20;
NX=100;
k1=[-30:2:-1];
k1=repmat(k1,1,10);
k2=-10;
dt=0.01;
eqstart1=1;


tic
% parpool(40)


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


a = -sqrt(3);  % 下界
b = sqrt(3);   % 上界
for i=1:10
Omegas(i,:) = (b - a) * rand(1, N) + a;
end


parfor i=1:cnt 
    disp(i)
    tt=ceil(i/15);
    Omega=Omegas(tt,:);
    [R1(i),RR1(i),RB(i)] = kurad2(tmax,kgd1,kgd2,edges,tri,...
        Omega,Ck1.Value(i),Ck2.Value(i),dt,eqstart1);   
end

save ome_uni10.mat

toc
