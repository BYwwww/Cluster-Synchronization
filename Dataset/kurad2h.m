function [RR1,RR2,theta_new] = kurad2h(tmax,kgd1,kgd2,edges,tri,Omega,K_1,K_2,dt,eqstart)
%% kuramoto
% tmax=2000;
% dt=0.01;
mms=fix(tmax/2/dt);
kmean1=mean(kgd1);
kmean2=mean(kgd2);
N=length(kgd2);
K1=K_1;

lenk=length(K_2);


% \theta0
if eqstart==1
    theta_now=zeros(N,1);  % now
else
    theta_now=rand(N,1)*2*pi;
end
theta_new=zeros(N,1);


% R = zeros(1,lenk);  %R：R-t series
% R(:,1)=abs(mean(exp(theta_now*1i)));

RR1=zeros(1,tmax/dt);
RR2=zeros(1,tmax/dt);
% R2 = zeros(1,tmax/dt);  %R2：R2-t series
% R2(1)=abs(mean(exp(theta_now*2i)));


for times=1:lenk
    
    K2=K_2(times);
%     R(1,times)=abs(mean(exp(theta_now*1i)));
    cnt=1;

    for t=0:dt:tmax  

        for n=1:N  %every node    
            theta_new(n) = theta_now(n) + (Omega(n) + K1 * delta_1(n)/kmean1 + K2 * delta_2(n)/kmean2)*dt; % 更新theta
        end

        % update R
        RR1(cnt)=abs(mean(exp(theta_new*1i)));       
        RR2(cnt)=abs(mean(exp(theta_new*2i)));
        cnt=cnt+1;
%         restnum=(tmax-t)/dt;
%         if(restnum<=mms-1)
%             RR1(times)=RR1(times)+abs(mean(exp(theta_new*1i)))/mms;
%             RR2(times)=RR2(times)+abs(mean(exp(theta_new*2i)))/mms;
%         end
        
        theta_now=theta_new; 
       
    end
    
    if mod(times,10)==0
        disp(['ten...']);

    end

end



function [d_theta2]=delta_2(node)
    len=kgd2(node);  %triangles
    if len==0
        d_theta2=0;
    else 
%         ids=squeeze([tri(node,1,:);tri(node,2,:)]);
        ids=reshape(tri(node,:,1:len),2,len);
        ids1=ids(1,:);  %node1
        ids2=ids(2,:);  %node2
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

