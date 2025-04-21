function [N,A,kgd1,kgd2,edges,tris,trin] = dataReal2()

A=load('hypertext2009_w.mat');
A=A.w;
A=full(A);
kgd1=sum(A);
N=size(A,1);

for i=1:N
    aa=find(A(i,:)==1);  %neighbors for subgraph       
    m=length(aa); 
    edges(i,1:m)=aa;
end


tris=load('hypertext2009_triangles.mat');
tris=tris.triangles;

trin=zeros(N,2);
kgd2=zeros(1,N);

for i=1:size(tris,1)
    tri_now=tris(i,:);
    i1=tri_now(1);
    i2=tri_now(2);
    i3=tri_now(3);
    Triangle2(i1,i2,i3);
end



function Triangle2(i1,i2,i3)
    kgd2([i1,i2,i3])=kgd2([i1,i2,i3])+1;  %triangles counted
    trin(i1,1,kgd2(i1))=i2;
    trin(i1,2,kgd2(i1))=i3;
    trin(i2,1,kgd2(i2))=i1;
    trin(i2,2,kgd2(i2))=i3;
    trin(i3,1,kgd2(i3))=i1;
    trin(i3,2,kgd2(i3))=i2;
end

end
