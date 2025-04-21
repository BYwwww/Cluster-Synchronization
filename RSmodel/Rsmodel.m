function [kg2,trin] = Rsmodel(N,A,kg20,tri1,alpha)
%% RSmodel
% Record a list of triangles for the network, delete some triangles, and randomly pick the same number of triangles 
% (equivalent to random coloring)

%------------------
% M:egdes 
% ntris: total number of triangles 
% ------------------ 
% A: adjacency matrix 
% kg: first order degree 
% kg2: generalized degree 
% edges: list of neighbors for each node 
% trin: list of triangles

kg2=kg20;
% A=sparse(N,N);
sumx=sum(kg20)/3;
tris=zeros(sumx,3);
cnt=1;

for j=1:N
    len=kg20(j);
    if len>0
        ids=reshape(tri1(j,:,1:len),2,len);
        ids1=ids(1,:);  %node1
        ids2=ids(2,:);  %node2
        for k=1:len
            id1=ids1(k);
            id2=ids2(k);
            tri_now=sort([j,id1,id2]);
            if ismember(tri_now,tris,'rows')==0    %不在原有三角形列表里
                tris(cnt,:)=tri_now;            %加入原有列表           
                cnt=cnt+1;
            end
        end
    end
end

% tris=unique(tris,'rows','sorted');
% cntt=size(tris,1)
cnt=cnt-1;
tris_change=fix(sumx*alpha);   %rewiring
ct=0;
while ct<tris_change
    xx=fix(cnt*rand()+1);
    tris_no=tris(xx,:);
    kg2(tris_no)=kg2(tris_no)-1;  %delete the label
    tris(xx,:)=[];
    cnt=cnt-1;
    ct=ct+1;
end

%Randomly generate triangles 
triss=zeros(sumx,3);
cntt=1;
%all triangles
for i=1:N
    aa=find(A(i,:)==1);  %Find neighbor nodes of subgraphs       
    m=length(aa); 
    if m>1
       B=triu(A(aa,aa));  % Extracting Upper Triangle Matrix of Subgraphs
       for j=1:m
           for k=j+1:m
               if B(j,k)==1
                   tri_now=sort([i,aa(j),aa(k)]);
                   triss(cntt,:)=tri_now;            %Add to the list of all triangles          
                   cntt=cntt+1;
               end
           end
       end
    end
end
% cntts=cntt-1;
trisss=unique(triss,'rows','sorted');
trittt=size(trisss,1);
[trins,~] = setdiff(triss,tris,'rows'); %De-duplication, de-existing
cntt=size(trins,1);  %Number of triangles left for selection
a=randperm(cntt);
aa=a(1:tris_change);  %Random selection of triangles
bb=trins(aa,:);
cc=[tris;bb];
trin=zeros(N,2);
kg2=zeros(1,N);

for i=1:sumx
    tri_now=cc(i,:);
    i1=tri_now(1);
    i2=tri_now(2);
    i3=tri_now(3);
    Triangle2(i1,i2,i3);
end


%***********************************************


function Triangle2(i1,i2,i3)
    kg2([i1,i2,i3])=kg2([i1,i2,i3])+1;  %现有三角形数
    trin(i1,1,kg2(i1))=i2;
    trin(i1,2,kg2(i1))=i3;
    trin(i2,1,kg2(i2))=i1;
    trin(i2,2,kg2(i2))=i3;
    trin(i3,1,kg2(i3))=i1;
    trin(i3,2,kg2(i3))=i2;
end


end