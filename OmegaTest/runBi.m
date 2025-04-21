clear;
N=1000;
tmax=1000;     %Ê±³¤
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


omega0 = 0.5;    % center
Delta = 1;       % std
n_samples = 1000; 
% Sampling from two normal distributions
samples = zeros(1, n_samples); % initialize
for j = 1:10
for i = 1:n_samples
    % Randomly select which peak (50% probability)
    if rand > 0.5
        samples(i) = normrnd(omega0, Delta);  % N(omega0, Delta^2) 
    else
        samples(i) = normrnd(-omega0, Delta); % N(-omega0, Delta^2) 
    end
end
Omegas(j,:)=samples;
end


parfor i=1:cnt  
    tt=ceil(i/15);
    Omega=Omegas(tt,:);
    [R1(i),RR1(i),RB(i)] = kurad2(tmax,kgd1,kgd2,edges,tri,...
        Omega,Ck1.Value(i),Ck2.Value(i),dt,eqstart1);   
end

save ome_Bi10.mat

toc
