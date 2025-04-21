function [RR1,RR2,theta_his] = kurad2_timess(tmax,kgd1,kgd2,edges,tri,Omega,K_1,K_2,dt,eqstart)
%% kuramoto，save history theta
% tmax=2000;
% dt=0.01;
mms=fix(tmax/3/dt);
kmean1=mean(kgd1);
kmean2=mean(kgd2);
N=length(kgd2);
K2=K_2;
K1=K_1;


%\theta0
%theta_now=eqstart;
if eqstart==1
    theta_now=zeros(N,1);  
else
    theta_now=rand(N,1)*2*pi;
end
theta_new=zeros(N,1);


at=fix(tmax/dt);
RR1=zeros(at,1);
RR2=zeros(at,1);
theta_his=zeros(N,1e4);

cnt=1;
cntt=1;
for t=0:dt:tmax
    for n=1:N  
        theta_new(n) = theta_now(n) + (Omega(n) + K1 * delta_1(n)/kmean1 + K2 * delta_2(n)/kmean2)*dt; % 更新theta
    end
    
    % update R
    RR1(cnt)=abs(mean(exp(theta_new*1i)));
    RR2(cnt)=abs(mean(exp(theta_new*2i)));
    if (tmax-t)/dt<1e5
    theta_his(:,cntt)=theta_new;
    cntt=cntt+1;
    end
    cnt=cnt+1;
   
    theta_now=theta_new;

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

