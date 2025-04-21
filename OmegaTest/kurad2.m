function [RR1,RR2,RB] = kurad2(tmax,kgd1,kgd2,edges,tri,Omega,K_1,K_2,dt,eqstart)
%% kuramoto
% tmax=2000;
% dt=0.01;
mms=fix(tmax/3/dt);
kmean1=mean(kgd1);
kmean2=mean(kgd2);
N=length(kgd2);
K1=K_1;

lenk=length(K_2);


%初相 \theta0
if eqstart==1
    theta_now=zeros(N,1);  
else
    theta_now=rand(N,1)*2*pi;
end
theta_new=zeros(N,1);


R = zeros(1,lenk);  %R：R-t
% R(:,1)=abs(mean(exp(theta_now*1i)));

RR1=zeros(1,lenk);
RR2=zeros(1,lenk);
RB=zeros(1,lenk);
% R2 = zeros(1,tmax/dt);  %R2：R2-t
% R2(1)=abs(mean(exp(theta_now*2i)));

% Finding subscripts for positive and negative omega
positive_indices = Omega > 0;  
negative_indices = Omega < 0;  


for times=1:lenk
    
    K2=K_2(times);
    cnt=1;

    for t=0:dt:tmax  

        for n=1:N  %every node    
            theta_new(n) = theta_now(n) + (Omega(n) + K1 * delta_1(n)/kmean1 + K2 * delta_2(n)/kmean2)*dt; % 更新theta
        end

        restnum=(tmax-t)/dt;
        if(restnum<=mms-1)
            RR1(times)=RR1(times)+abs(mean(exp(theta_new*1i)))/mms;
            RR2(times)=RR2(times)+abs(mean(exp(theta_new*2i)))/mms;
            RBtemp=abs(mean(exp(theta_new(positive_indices)*1i)))/2+abs(mean(exp(theta_new(negative_indices)*1i)))/2-abs(mean(exp(theta_new*1i)));
            RB(times)=RB(times)+RBtemp/mms;
        end
        
        theta_now=theta_new; 
       
    end
    
    if mod(times,10)==0
        disp(['ten...']);
    end

end



function [d_theta2]=delta_2(node)
    len=kgd2(node);  
    if len==0
        d_theta2=0;
    else 
%         ids=squeeze([tri(node,1,:);tri(node,2,:)]);
        ids=reshape(tri(node,:,1:len),2,len);
        ids1=ids(1,:);  
        ids2=ids(2,:);  
        d_theta2=sum(sin(theta_now(ids1)+theta_now(ids2)-2*theta_now(node)));  % Σsin(j+k-2i)
    end
end


function [d_theta1]=delta_1(node)   
    len1=kgd1(node);
    if len1==0
        d_theta1=0;
    else
        ids=edges(node,1:len1);
        d_theta1=sum(sin(theta_now(ids)-theta_now(node))); % Σsin(j-i)
    end
end





end

