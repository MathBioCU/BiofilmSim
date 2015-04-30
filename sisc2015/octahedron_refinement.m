function [num_t,num_v,num_e,v,t_to_v,v_to_t,e_to_v,v_to_e,v_to_v] = octahedron_refinement( ref)
%This program starts off with an octahedron and calculates a number of
%different matrices that relate to the connectivity and position of edges,
%vertices, and faces. Then, the octahedron and all of the data calculated
%with it is refined to include more and more points.
%
%A note on notation:  X{position in cell array}(row)(column) 
%is commonly used in this program. The position in the cell array notes on 
%which refinement level the program is working on. It is necessary not to 
%rewrite the data for lower refinement levels than the current since that
%data is %necessary in the multigrid algorithm later on in the main 
%program.
%**************************************************************************

%The first section of this program is simply done by hand.
%By definition an octahedron has 8 triangular faces
%(num_t=number of triangles), 6 vertices (num_v=number of vertices), and 12
%edges (num_e=number of edges). num_t, num_v, and num_e are vectors of
%length ref+1 because these values are stored for each level of
%refinement. Because matlab numbers vectors and matrices starting from 1,
%they must be of length ref+1. The first value does not correspond
%with the first refinement but the octahedron, unrefined.
num_t=zeros(ref+1,1); num_v=num_t; num_e=num_t; 

num_t(1)=8; num_v(1)=6; num_e(1)=12;

%these next cell arrays are set up to store different data from each
%refinement level. Each entry in the cell will be an array of varying size
%corresponding to the octahedron and each subsequent refinement.
t_to_v=cell(ref+1,1); 
v_to_t=t_to_v; e_to_t=t_to_v; t_to_e=t_to_v;
e_to_v=t_to_v; v_to_e=t_to_v; v_to_v=t_to_v;

%Because of the way vertices will be numbered, it is not necessary to
%rewrite old vertices and change the numbering for each refinement.
%Instead, the first 6 vertices will always be the first 6 for each
%refinement etc.

%vertices
v = [  1  0  0 ;
       0  1  0 ;
      -1  0  0 ;        
       0 -1  0 ;
       0  0  1 ;
       0  0 -1 ];

v_to_t{1}=zeros(6,6);
ipt=ones(6,1);

%triangle to vertex and vertex to triangle
n=1; p=1; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=1; p=2; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=1; p=5; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=2; p=2; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=2; p=3; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=2; p=5; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=3; p=3; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=3; p=4; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=3; p=5; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=4; p=4; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=4; p=1; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=4; p=5; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=5; p=2; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=5; p=1; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=5; p=6; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=6; p=3; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=6; p=2; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=6; p=6; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=7; p=4; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=7; p=3; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=7; p=6; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=8; p=1; t_to_v{1}(n,1)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=8; p=4; t_to_v{1}(n,2)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=8; p=6; t_to_v{1}(n,3)=p; v_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;

%edge to vertex and vertex to edge

e_to_v{1}=zeros(12,2);
v_to_e{1}=zeros(6,6);
ipt=ones(6,1);

n=1; p=1; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=1; p=5; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=2; p=2; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=2; p=5; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=3; p=3; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=3; p=5; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=4; p=4; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=4; p=5; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=5; p=1; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=5; p=6; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=6; p=2; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=6; p=6; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=7; p=3; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=7; p=6; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=8; p=4; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=8; p=6; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=9; p=1; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=9; p=2; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=10; p=2; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=10; p=3; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=11; p=3; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=11; p=4; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=12; p=4; e_to_v{1}(n,1)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 
n=12; p=1; e_to_v{1}(n,2)=p; v_to_e{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1; 

%triangle to edge and edge to triangle
t_to_e{1}=zeros(8,3);
e_to_t{1}=zeros(12,2);
ipt=ones(12,1);

n=1; p=9; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=1; p=2; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=1; p=1; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=2; p=10; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=2; p=3; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=2; p=2; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=3; p=11; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=3; p=4; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=3; p=3; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=4; p=12; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=4; p=1; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=4; p=4; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=5; p=9; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=5; p=5; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=5; p=6; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=6; p=10; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=6; p=6; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=6; p=7; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=7; p=11; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=7; p=7; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=7; p=8; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=8; p=12; t_to_e{1}(n,1)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=8; p=8; t_to_e{1}(n,2)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;
n=8; p=5; t_to_e{1}(n,3)=p; e_to_t{1}(p,ipt(p))=n; ipt(p)=ipt(p)+1;

%adjacent vertices
v_to_v{1}=zeros(6,6);
v_to_v{1}(1,1)=2; v_to_v{1}(1,2)=5; v_to_v{1}(1,3)=4; v_to_v{1}(1,4)=6;
v_to_v{1}(2,1)=3; v_to_v{1}(2,2)=5; v_to_v{1}(2,3)=1; v_to_v{1}(2,4)=6;
v_to_v{1}(3,1)=4; v_to_v{1}(3,2)=5; v_to_v{1}(3,3)=2; v_to_v{1}(3,4)=6;
v_to_v{1}(4,1)=1; v_to_v{1}(4,2)=5; v_to_v{1}(4,3)=3; v_to_v{1}(4,4)=6;
v_to_v{1}(5,1)=1; v_to_v{1}(5,2)=2; v_to_v{1}(5,3)=3; v_to_v{1}(5,4)=4;
v_to_v{1}(6,1)=1; v_to_v{1}(6,2)=2; v_to_v{1}(6,3)=3; v_to_v{1}(6,4)=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%refinement

for ijk=1:ref
    
    L=num_e(ijk);
    
    %these formulas for the number of triangles, edges, and vertices can be
    %proven. They all increase at a rate of approximately n_new=4*n_old.
    t_to_e_new=zeros(4*num_t(ijk),3); v_to_e{ijk+1}=zeros(num_v(ijk)+L,6); t_to_v_new=zeros(num_t(ijk)*4,3);
    v_to_t{ijk+1}=zeros(num_v(ijk)+L,6); e_to_v_new=zeros(4*num_e(ijk),2); e_to_t{ijk+1}=zeros(num_e(ijk)*4,2);
    
    %these are pointers that will be used later to mark the left-most column of some
    %the arrays that is empty
    ipt=ones(4*num_t(ijk),1); ipt2=ones(4*num_e(ijk),1); ipt3=ones(4*num_e(ijk),1); ipt4=ones(num_v(ijk)+L,1);
    
    v_new=zeros(num_v(ijk)+L,3);
    %All of the vertices from previous refinements remain unchanged in
    %location and pointer number
    v_new(1:num_v(ijk),:)=v;
    %New vertices are formed by cycling through the edges and forming a new
    %vertex at the center of each edge.
    for i=1:num_e(ijk),
        v1=v(e_to_v{ijk}(i,1),:);
        v2=v(e_to_v{ijk}(i,2),:);
        v_new(num_v(ijk)+i,:)=(v1+v2)/2;
    end;
    v=v_new;  
    
    %edge to vertex and vertex to edge (part of it)
    
    %each of the old edges now has a vertex in the center of it. So, there
    %are now twice as many edges. The new edges are labaled such that the
    %half of the old edge nearer the first old vertex is now edge 2*i-1 and
    %the other half of the old edge is now edge 2*i
    %
    %the vertex to edge to data is filled simultaneously but not in any
    %clear ordering. The ipt2 vector is used to prevent rewriting of data
    %in this array. If e_to_v(i,?)=k, v_to_e(k,?)=i, they are inverses.
    for i=1:num_e(ijk)
        e_to_v_new(2*i-1,1)=e_to_v{ijk}(i,1); v_to_e{ijk+1}(e_to_v{ijk}(i,1),ipt2(e_to_v{ijk}(i,1)))=2*i-1; ipt2(e_to_v{ijk}(i,1))=ipt2(e_to_v{ijk}(i,1))+1;
        e_to_v_new(2*i-1,2)=i+num_v(ijk);     v_to_e{ijk+1}(i+num_v(ijk),ipt2(i+num_v(ijk)))=2*i-1;         ipt2(i+num_v(ijk))=ipt2(i+num_v(ijk))+1;

        e_to_v_new(2*i,1)=i+num_v(ijk);       v_to_e{ijk+1}(i+num_v(ijk),ipt2(i+num_v(ijk)))=2*i;           ipt2(i+num_v(ijk))=ipt2(i+num_v(ijk))+1;
        e_to_v_new(2*i,2)=e_to_v{ijk}(i,2);   v_to_e{ijk+1}(e_to_v{ijk}(i,2),ipt2(e_to_v{ijk}(i,2)))=2*i;   ipt2(e_to_v{ijk}(i,2))=ipt2(e_to_v{ijk}(i,2))+1;        
    end
    
    for i=1:num_t(ijk)
        
        i1=t_to_v{ijk}(i,1);  j1=t_to_e{ijk}(i,1);   k1=j1+num_v(ijk);
        i2=t_to_v{ijk}(i,2);  j2=t_to_e{ijk}(i,2);   k2=j2+num_v(ijk);
        i3=t_to_v{ijk}(i,3);  j3=t_to_e{ijk}(i,3);   k3=j3+num_v(ijk); 
        
         %fill in simultaneously
         
        %each old triangle will be broken up into 4 new ones when the edges
        %are cut in half. This section just represents a lot of book
        %keeping to make sure that the write triangles match up with the
        %right vertices. Triangle vetices are labeled counter clockwise.
        t_to_v_new((i-1)*4+1,1)=i1; v_to_t{ijk+1}(i1,ipt(i1))=4*(i-1)+1; ipt(i1)=ipt(i1)+1;
        t_to_v_new((i-1)*4+2,1)=i2; v_to_t{ijk+1}(i2,ipt(i2))=4*(i-1)+2; ipt(i2)=ipt(i2)+1;
        t_to_v_new((i-1)*4+3,1)=i3; v_to_t{ijk+1}(i3,ipt(i3))=4*(i-1)+3; ipt(i3)=ipt(i3)+1;
        
        t_to_v_new(4*i,1)=k1; v_to_t{ijk+1}(k1,ipt(k1))=4*i; ipt(k1)=ipt(k1)+1;
        t_to_v_new(4*i,2)=k2; v_to_t{ijk+1}(k2,ipt(k2))=4*i; ipt(k2)=ipt(k2)+1;
        t_to_v_new(4*i,3)=k3; v_to_t{ijk+1}(k3,ipt(k3))=4*i; ipt(k3)=ipt(k3)+1;
                          
        t_to_v_new((i-1)*4+1,2)=k1; v_to_t{ijk+1}(k1,ipt(k1))=4*(i-1)+1; ipt(k1)=ipt(k1)+1;
        t_to_v_new((i-1)*4+1,3)=k3; v_to_t{ijk+1}(k3,ipt(k3))=4*(i-1)+1; ipt(k3)=ipt(k3)+1;
        t_to_v_new((i-1)*4+2,2)=k2; v_to_t{ijk+1}(k2,ipt(k2))=4*(i-1)+2; ipt(k2)=ipt(k2)+1;
        t_to_v_new((i-1)*4+2,3)=k1; v_to_t{ijk+1}(k1,ipt(k1))=4*(i-1)+2; ipt(k1)=ipt(k1)+1;
        t_to_v_new((i-1)*4+3,2)=k3; v_to_t{ijk+1}(k3,ipt(k3))=4*(i-1)+3; ipt(k3)=ipt(k3)+1;
        t_to_v_new((i-1)*4+3,3)=k2; v_to_t{ijk+1}(k2,ipt(k2))=4*(i-1)+3; ipt(k2)=ipt(k2)+1;
        
        %the remaining edge to vertex and vertex to edge data is now filled. This data all
        %comes from new edges connecting two new vertices.
        e_to_v_new(3*(i-1)+1+num_e(ijk)*2,1)=k1; v_to_e{ijk+1}(k1,ipt2(k1))=3*(i-1)+1+num_e(ijk)*2; ipt2(k1)=ipt2(k1)+1;
        e_to_v_new(3*(i-1)+1+num_e(ijk)*2,2)=k2; v_to_e{ijk+1}(k2,ipt2(k2))=3*(i-1)+1+num_e(ijk)*2; ipt2(k2)=ipt2(k2)+1;

        e_to_v_new(3*(i-1)+2+num_e(ijk)*2,1)=k2; v_to_e{ijk+1}(k2,ipt2(k2))=3*(i-1)+2+num_e(ijk)*2; ipt2(k2)=ipt2(k2)+1;
        e_to_v_new(3*(i-1)+2+num_e(ijk)*2,2)=k3; v_to_e{ijk+1}(k3,ipt2(k3))=3*(i-1)+2+num_e(ijk)*2; ipt2(k3)=ipt2(k3)+1;

        e_to_v_new(3*(i-1)+3+num_e(ijk)*2,1)=k3; v_to_e{ijk+1}(k3,ipt2(k3))=3*(i-1)+3+num_e(ijk)*2; ipt2(k3)=ipt2(k3)+1;
        e_to_v_new(3*(i-1)+3+num_e(ijk)*2,2)=k1; v_to_e{ijk+1}(k1,ipt2(k1))=3*(i-1)+3+num_e(ijk)*2; ipt2(k1)=ipt2(k1)+1;
       
        % now the edge to triangle and triangle to edge data can be determined. This is the most complex part, 
        % but also the part that is used to least later on in the program.
        % The complexity comes from the fact that each triangle can end up
        % with sort of a different orientation in some sense and it is
        % necessary to determine what this orientation is to known what
        % edges and what triangles are touching.
        for j=1:2
            e1=v_to_e{ijk+1}(k1,j);
            e2=v_to_e{ijk+1}(k2,j);
            e3=v_to_e{ijk+1}(k3,j);
            if e_to_v_new(e1,1)==i1 || e_to_v_new(e1,2)==i1
                t_to_e_new(4*(i-1)+1,1)=e1;  e_to_t{ijk+1}((e1),ipt3(e1))=4*(i-1)+1; ipt3(e1)=ipt3(e1)+1;
            end
            if e_to_v_new(e3,1)==i1 || e_to_v_new(e3,2)==i1
                t_to_e_new(4*(i-1)+1,3)=e3;  e_to_t{ijk+1}(e3,ipt3(e3))=4*(i-1)+1; ipt3(e3)=ipt3(e3)+1;
            end
            if e_to_v_new(e1,1)==i2 || e_to_v_new(e1,2)==i2
                t_to_e_new(4*(i-1)+2,3)=e1;  e_to_t{ijk+1}(e1,ipt3(e1))=4*(i-1)+2; ipt3(e1)=ipt3(e1)+1;
            end
            if e_to_v_new(e2,1)==i2 || e_to_v_new(e2,2)==i2
                t_to_e_new(4*(i-1)+2,1)=e2;  e_to_t{ijk+1}(e2,ipt3(e2))=4*(i-1)+2; ipt3(e2)=ipt3(e2)+1;
            end
            if e_to_v_new(e2,1)==i3 || e_to_v_new(e2,2)==i3
                t_to_e_new(4*(i-1)+3,3)=e2;  e_to_t{ijk+1}(e2,ipt3(e2))=4*(i-1)+3; ipt3(e2)=ipt3(e2)+1;
            end
            if e_to_v_new(e3,1)==i3 || e_to_v_new(e3,2)==i3
                t_to_e_new(4*(i-1)+3,1)=e3;  e_to_t{ijk+1}(e3,ipt3(e3))=4*(i-1)+3; ipt3(e3)=ipt3(e3)+1;
            end
        end
        t_to_e_new(4*(i-1)+1,2)=3*(i-1)+3+num_e(ijk)*2; e_to_t{ijk+1}(3*(i-1)+3+num_e(ijk)*2,ipt3(3*(i-1)+3+num_e(ijk)*2))=4*(i-1)+1; ipt3(3*(i-1)+3+num_e(ijk)*2)=ipt3(3*(i-1)+3+num_e(ijk)*2)+1;       
        t_to_e_new(4*(i-1)+2,2)=3*(i-1)+1+num_e(ijk)*2; e_to_t{ijk+1}(3*(i-1)+1+num_e(ijk)*2,ipt3(3*(i-1)+1+num_e(ijk)*2))=4*(i-1)+2; ipt3(3*(i-1)+1+num_e(ijk)*2)=ipt3(3*(i-1)+1+num_e(ijk)*2)+1;
        t_to_e_new(4*(i-1)+3,2)=3*(i-1)+2+num_e(ijk)*2; e_to_t{ijk+1}(3*(i-1)+2+num_e(ijk)*2,ipt3(3*(i-1)+2+num_e(ijk)*2))=4*(i-1)+3; ipt3(3*(i-1)+2+num_e(ijk)*2)=ipt3(3*(i-1)+2+num_e(ijk)*2)+1;
       
        t_to_e_new(4*(i-1)+4,1)=3*(i-1)+1+num_e(ijk)*2; e_to_t{ijk+1}(3*(i-1)+1+num_e(ijk)*2,ipt3(3*(i-1)+1+num_e(ijk)*2))=4*(i-1)+2; ipt3(3*(i-1)+1+num_e(ijk)*2)=ipt3(3*(i-1)+1+num_e(ijk)*2)+1;
        t_to_e_new(4*(i-1)+4,2)=3*(i-1)+2+num_e(ijk)*2; e_to_t{ijk+1}(3*(i-1)+2+num_e(ijk)*2,ipt3(3*(i-1)+2+num_e(ijk)*2))=4*(i-1)+3; ipt3(3*(i-1)+2+num_e(ijk)*2)=ipt3(3*(i-1)+2+num_e(ijk)*2)+1;
        t_to_e_new(4*(i-1)+4,3)=3*(i-1)+3+num_e(ijk)*2; e_to_t{ijk+1}(3*(i-1)+3+num_e(ijk)*2,ipt3(3*(i-1)+3+num_e(ijk)*2))=4*(i-1)+1; ipt3(3*(i-1)+3+num_e(ijk)*2)=ipt3(3*(i-1)+3+num_e(ijk)*2)+1;    
        
    end
    
    %find connected vertices
    %this data can be determined just from the previous edge to
    %vertex data already known from before.
    v_to_v_new=zeros(num_v(ijk+1),6);
    for i=1:num_e(ijk)*4
        v1=e_to_v_new(i,1);
        v2=e_to_v_new(i,2);
        v_to_v_new(v1,ipt4(v1))=v2; ipt4(v1)=ipt4(v1)+1; 
        v_to_v_new(v2,ipt4(v2))=v1; ipt4(v2)=ipt4(v2)+1; 
    end
    e_to_v{ijk+1}=e_to_v_new;
    t_to_v{ijk+1}=t_to_v_new;
    t_to_e{ijk+1}=t_to_e_new;
    v_to_v{ijk+1}=v_to_v_new;
    
    %determine the number of edges, triangles and vertices for the next
    %refinement.
    num_v(ijk+1)=num_v(ijk)+num_e(ijk);
    num_t(ijk+1)=num_t(ijk)*4; 
    num_e(ijk+1)=num_e(ijk)*4;
 
    %project
    %make each vertex a vector of length 1, this guarantees that it will be
    %on the surface of the unit sphere.
    for i=1:length(v(:,1))
        v_new(i,1)=v(i,1)/norm(v(i,:));
        v_new(i,2)=v(i,2)/norm(v(i,:));
        v_new(i,3)=v(i,3)/norm(v(i,:));
    end
    v=v_new;
    
end

end

