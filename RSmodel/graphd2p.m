function [A,kgd1,kgd2,edges,tri] = graphd2p(N,m,NX)
%% Generate pure second-order simplicial complexes

Avoid=1;  %allow rollback
maxn=(N-1)*(N-2)*0.5;
A=sparse(N,N);
% generalized degree
kgi=zeros(1,N);


%main   
xaus=4;
while xaus>3
    %initialize 
    disp(['Init...']);
%         for i=1:N
%             kgi(i)=fix(m*(rand()^(-1/(gamma-1))));  %power law           
%             while kgi(i)>maxn
%                 kgi(i)=fix(m*(rand()^(-1/(gamma-1))));
%             end
%         end
    kgi=fix(poissrnd(m,[1,N]));
%     kgi=random('Poisson',m,[1,N]);
%     kgi=fix(randn(1,N)*m/6+m);
    
    kkk=max(kgi);
    kgd2=zeros(1,N);
    A=sparse(N,N);
    tri=zeros(N,2,kkk);

    xaus=sum(kgi);  %all threads
    naus=0; % max illegal times 

    while (xaus>3)&&(naus<1+Avoid*NX)
        %1
        x1=xaus*rand();
        i1=Choose(x1);
        kgi(i1)=kgi(i1)-1; % unpaired threads-1
        xaus=xaus-1;

        %2
        x2=xaus*rand();
        i2=Choose(x2);
        kgi(i2)=kgi(i2)-1;
        xaus=xaus-1;

        %3
        x3=xaus*rand();
        i3=Choose(x3);
        kgi(i3)=kgi(i3)-1;
        xaus=xaus-1;

        kgd2([i1,i2,i3])=kgd2([i1,i2,i3])+1;  % paired threads+1

        %Check if a triangle already exists
        if (i1~=i2)&&(i1~=i3)&&(i2~=i3)&&(Check(i1,i2,i3)==0)
            %legal，Generating Triangle
            Triangle(i1,i2,i3);
            A(i1,i2)=1;
            A(i2,i1)=1;
            A(i1,i3)=1;
            A(i3,i1)=1;
            A(i3,i2)=1;
            A(i2,i3)=1;
        else
            %illegal，back+1
            naus=naus+1;
            if Avoid==1
                kgd2([i1,i2,i3])=kgd2([i1,i2,i3])-1;
                kgi([i1,i2,i3])=kgi([i1,i2,i3])+1;
%                     xaus=xaus+3;
            end
        end      
    end
end



%one order information
kgd1=zeros(1,N);
for j=1:N
    kgd1(j)=sum(A(j,:));
end
edges=zeros(N,max(kgd1));

for j=1:N
    ee=find(A(j,:)==1);
    edges(j,1:kgd1(j))=ee;
end


% plot(graph(A));

%*************************************************************************
%random selection
function [node]=Choose(x)
    cdfx=cumsum(kgi);
    node=find(x<cdfx,1);
end

%Check if a triangle already exists
function [flag]=Check(i1,i2,i3)
    flag=0;
    for j=1:kgd2(i1)
        if (tri(i1,1,j)==i2&&tri(i1,2,j)==i3)||(tri(i1,1,j)==i3&&tri(i1,2,j)==i2)
            flag=1;
            break;
        end
    end
end

%Generating Triangles
function Triangle(i1,i2,i3)
    tri(i1,1,kgd2(i1))=i2;
    tri(i1,2,kgd2(i1))=i3;
    tri(i2,1,kgd2(i2))=i1;
    tri(i2,2,kgd2(i2))=i3;
    tri(i3,1,kgd2(i3))=i1;
    tri(i3,2,kgd2(i3))=i2;
end

end

