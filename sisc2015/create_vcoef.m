% close all
% clear all
% clc

% numtimesteps=25;

levelsV=5;
levelsP=4;

%%%%%%%%%%%%%dimensional constants
% w=1;
e0=4*1.8*10*10^-6*tan(0.13);
v0=e0*w;%speed at the middle
visc0=10^-3;%this is dynamic viscosity of water, units of kg/m/s
rho0=998;%9.983*10^(-7);%in units of kg/m^3
charLength=10*10^-6;
timefreq=w;
d0mean=1.5921*10^-6;%calculated from given Stewart data as (Vol/(# Lagr nodes interior to the Vol))^(1/3) in the M-file Determine_d0.m
d0mean=d0mean/charLength;%nondimensionalized to use in the simulations

%%%%%%%%%%%%%%%%%%%%%%

mult=10^-6;

% xlength=10*mult;%diameter of the tube
% ylength=10*mult;%radius of the tube
% zlength=10*mult;%length of the tube for the computational domain

% Vmax=laminarVel(xlength/2,ylength/2,xlength,ylength);
% pgrad=-v0*visc0/xlength^2/Vmax;
% LaminarInit=interpolate2htoh(LaminarInit);
% ux(:,1)=LaminarInit;
% p0=-pgrad*zlength;
p0=1;

%%%%%%%%%%%%nondimensional constants
st=timefreq*charLength/v0;
eu=p0/rho0/v0^2;
Re=rho0*charLength*v0/visc0;
% fmax=1;
% f0=175000;%3 is because 3 times rest length is the breaking criteria
fc=charLength/rho0/v0^2;

eust=eu/st;


%%%%%%%%%%%%%%%%%%%%%%%%

% numtimesteps=100;

initdensity=1;
addldens=0.2;%
visc=1;
% addlvisc=500;

h=0.9/64;%64 should be min, this is the hx hy, and hz
% n=1/h-1;
% levels=5;
% dt=.01;
xlength=0.9;%diameter of the tube
ylength=1.8;%radius of the tube
zlength=0.9;%length of the tube for the computational domain
xvec=0:h:xlength;
yvec=0:h:ylength;
zvec=0:h:zlength;


[x,y,z]=meshgrid(xvec,yvec,zvec);

tend=numtimesteps*dt;
% connectdist=3/18;

% % load Stewart1.csv
% % load Stewart2.csv
% % % % % % % load Stewart3.csv
% % % % % % % Xtemp=[Stewart3(:,4),Stewart3(:,6),Stewart3(:,5)]*.06;
% % % % % % % 
% % % % % % % x1=Xtemp(Xtemp(:,1)>10&Xtemp(:,1)<22&Xtemp(:,2)>0&Xtemp(:,2)<16&Xtemp(:,3)>10&Xtemp(:,3)<23,:);
% % % % % % % % x2=Xtemp(Xtemp(:,1)>10&Xtemp(:,1)<19&Xtemp(:,2)>2&Xtemp(:,2)<8&Xtemp(:,3)>10&Xtemp(:,3)<19,:);
% % % % % % % 
% % % % % % % x1(:,1)=x1(:,1)-min(x1(:,1));
% % % % % % % x1(:,2)=x1(:,2)-min(x1(:,2));
% % % % % % % x1(:,3)=x1(:,3)-87/128*mean(x1(:,3));
% % % % % % % 
% % % % % % % % x4=lhsdesign(40,3);%25
% % % % % % % 
% % % % % % % % x1(:,1)=x1(:,1)*mult+10.5*mult;
% % % % % % % % x1(:,2)=x1(:,2)*mult-1.25*mult;
% % % % % % % % x1(:,3)=x1(:,3)*mult+60.5*mult;
% % % % % % % % x2(:,1)=x2(:,1)*mult+10.5*mult;
% % % % % % % % x2(:,2)=x2(:,2)*mult-2*mult;
% % % % % % % % x2(:,3)=x2(:,3)*mult+60.5*mult;
% % % % % % % % x3(:,1)=x3(:,1)*3*mult+23.5*mult;
% % % % % % % % x3(:,2)=x3(:,2)*3*mult+3*mult;
% % % % % % % % x3(:,3)=x3(:,3)*3*mult+73.5*mult;
% % % % % % % % x4(:,1)=x4(:,1)*9*mult+20.5*mult;
% % % % % % % % x4(:,3)=x4(:,3)*9*mult+70.5*mult;
% % % % % % % % x4(:,2)=0;
% % % % % % % % X=[x1k;x2k;x4k]/charLength;
%fID='/home/jstotsky/scratch/37_Deg_Test_1_Live_Cells.txt';
[Xb,Yb,Zb]=importfile('37 Deg Test 1, Live Cells.txt');
x1=[Xb, Yb,Zb];
x1=x1(x1(:,1)<20 & x1(:,2)<30 & x1(:,3)<15 & x1(:,2)>9 & x1(:,1)>10 ,:);
x1=[x1(:,1)-10,x1(:,2)-10,x1(:,3)]; %set edge from 10 microns to 0 microns
X=[x1(:,1),x1(:,3),x1(:,2)]*mult/charLength;

% X=X(find(X(:,2)<0.2),:);
X=X(X(:,1)<xlength,:);
X=X(X(:,2)>0,:);
% X=X(X(:,2)>0.01,:);
Xupper=X(X(:,3)>=ylength,:);
Xupper(:,3)=ylength;
X=X(X(:,3)<ylength,:);
X=[X;Xupper];

Xlower=X(X(:,3)<=0,:);
Xlower(:,3)=0;
X=X(X(:,3)>=0,:);
X=[X;Xlower];


% X=[0,0,0];


lower=find(X(:,3)==0);
upper=find(X(:,3)==ylength);

% load X.mat
% load X1.mat
[X1,Y1,Z1]=meshgrid(0:d0mean:0.2,0:d0mean:0.2,0:d0mean:0.2);
Xt=[X1(:),Y1(:),Z1(:)];
Xt=Xt*charLength;

X=[X(:,1),X(:,3),X(:,2)];

% X=[0.5 0.5 0.95]; Xupper=[0,0,0]; Xlower=[0,0,0];
% plot3(X(:,1),X(:,3),X(:,2),'r.')%,xvec,135*mult*ones(length(xvec)),[fliplr(yvec),yvec(2:length(yvec))])
% xlabel('x');
% ylabel('z');
% zlabel('y');
% hold on
% plot3(X1(:,1),X1(:,3),X1(:,2),'k.','MarkerSize',20)%,xvec,135*mult*ones(length(xvec)),[fliplr(yvec),yvec(2:length(yvec))])


% view(2)


numOfPoints=size(X,1);
numOfnonzero=numOfPoints-length(find(X(:,3)==0))-length(find(X(:,3)==1));
%% initialize all the variables
%A is connectivity matrix for Lagrangian points
A=zeros(numOfPoints);
d0=A;
% D=zeros(numOfPoints);
Xdist=A;
Ydist=A;
Zdist=A;
U=zeros(numOfPoints,3);
%next the initials for the eulerian points
XrelVel=Xdist;
YrelVel=Xdist;
ZrelVel=Xdist;
[Em,En,Ep]=size(x);



    
ux=zeros(Em,En,Ep);
uy=ux;
uz=ux;

% kw=(w/(2*visc0/rho0))^(1/2);
% for i=1:Em
%     uz(i,:,:)=sqrt(real(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2+imag(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2).*sin(angle(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))));
% end


p0=0;
p=zeros(Em,En,Ep);%constant pressure
% % p(:,:,1)=p0;
% % % p(Em,:,:)=pgrad*z(Em,:,:)+p0;
% % for c6=1:Em-1
% %     p(c6,:,:)=p(Em,:,:);
% % end

Efx=zeros(Em,En,Ep);
Efy=Efx;
Efz=Efx;


%% Set up the connections between the Lagrangian nodes (Biofilm nodes)
X(:,3)=mod(X(:,3),zlength);
X(:,1)=mod(X(:,1),xlength);
for c1=1:numOfPoints
    for c2=c1+1:numOfPoints
        dtemp=norm(X(c1,:)-X(c2,:));
        if dtemp<connectdist && dtemp>h/2
            
            A(c1,c2)=1;
            A(c2,c1)=1;

            d0(c1,c2)=dtemp;
            d0(c2,c1)=d0(c1,c2);
        end
    end
end
Xr=[X(:,1), X(:,2), X(:,3)+zlength];
for c1=1:numOfPoints
    for c2=c1+1:numOfPoints
        dtemp=norm(X(c1,:)-Xr(c2,:));
        if dtemp<connectdist && dtemp>h/2
            
            A(c1,c2)=1;
            A(c2,c1)=1;

            d0(c1,c2)=dtemp;
            d0(c2,c1)=d0(c1,c2);
        end
    end
end
Xl=[X(:,1), X(:,2), X(:,3)-zlength];
for c1=1:numOfPoints
    for c2=c1+1:numOfPoints
        dtemp=norm(X(c1,:)-Xl(c2,:));
        if dtemp<connectdist && dtemp>h/2
            
            A(c1,c2)=1;
            A(c2,c1)=1;

            d0(c1,c2)=dtemp;
            d0(c2,c1)=d0(c1,c2);
        end
    end
end

Xr=[X(:,1)+xlength, X(:,2), X(:,3)];
for c1=1:numOfPoints
    for c2=c1+1:numOfPoints
        dtemp=norm(X(c1,:)-Xr(c2,:));
        if dtemp<connectdist && dtemp>h/2
            
            A(c1,c2)=1;
            A(c2,c1)=1;

            d0(c1,c2)=dtemp;
            d0(c2,c1)=d0(c1,c2);
        end
    end
end
Xl=[X(:,1)-xlength, X(:,2), X(:,3)];
for c1=1:numOfPoints
    for c2=c1+1:numOfPoints
        dtemp=norm(X(c1,:)-Xl(c2,:));
        if dtemp<connectdist && dtemp>h/2
            
            A(c1,c2)=1;
            A(c2,c1)=1;

            d0(c1,c2)=dtemp;
            d0(c2,c1)=d0(c1,c2);
        end
    end
end

% gplotdc3D(A,50*X)


% axis(50*[0 xlength 0 ylength 0 zlength])
%     axis(chlng*[.35 .65 1.2 2.0 0 .2])
%    
   
%     xlabel('x')
%     ylabel('z')
%     zlabel('y')
%         view(50, 10) %side view
% b=zeros(numOfPoints,1);
% b=5;
b=b/e0./d0;
b(~A)=0;
t=0:dt:tend;
[rowind,colind]=find(A);
matind=find(A);

% m=densLagr*ones(numOfPoints,1);
% gplotdc(A,X)
% axis([-.5 w -.5 w])
% pause
Xdist;
Ydist;
% k=1*.042^2;%d0mean^2;%from the klapper paper, change to d0^2 for 3d
K=fmax./d0;
K(~A)=0;
% fcd0mean=fc/d0mean^3*K;
Xstore=X;
X0=X;
%Transfer the density to the Eulerian grid from the Lagrangian/compute
%new densities

Edens{1}.dx=h;
Edens{1}.Edensin=transferLtoEdens3Dper_e2(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);
viscmat{1}=transferLtoEvisc3Dper_e2(h,x,y,z,X,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        Edens{1}.Edensmidlr=(Edens{1}.Edensin(:,2:En,:)+Edens{1}.Edensin(:,1:En-1,:))/2;
        Edens{1}.Edensmidud=(Edens{1}.Edensin(2:Em,:,:)+Edens{1}.Edensin(1:Em-1,:,:))/2;
        Edens{1}.Edensmidfb=(Edens{1}.Edensin(:,:,2:Ep)+Edens{1}.Edensin(:,:,1:Ep-1))/2;
        
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
Edens{1}.x2h=x;
Edens{1}.y2h=y;
Edens{1}.z2h=z;

Edens{1}.Em=Em;
Edens{1}.En=En;
Edens{1}.Ep=Ep;
% save datatemp1.mat h hy x y z X Ind initdensity addldens d0mean

for c5=2:levelsV
    clear x2h y2h z2h dx hytemp yvectemp twotopow Em1 En1 Ep1 redf blackf red black
    twotopow=2^(c5-1);
    dx=h*twotopow;
    
    Edens{c5}.dx=dx;
    [x2h,y2h,z2h]=meshgrid(0:dx:xlength,0:dx:ylength,0:dx:zlength);
    
    [Em1,En1,Ep1]=size(x2h);
    
        
    Edens{c5}.x2h=x2h;
    Edens{c5}.y2h=y2h;
    Edens{c5}.z2h=z2h;
    Edens{c5}.Em=Em1;
    Edens{c5}.En=En1;
    Edens{c5}.Ep=Ep1;
   
% Edens{c5}.Edensin=transferLtoEdens3Dper(dx,x2h,y2h,z2h,X,initdensity,addldens,d0mean,numOfnonzero,zlength);
% viscmat{c5}=transferLtoEvisc3Dper(dx,x2h,y2h,z2h,X,visc,addlvisc,numOfnonzero,c5,zlength);
%         Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
%         Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
%         Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;
%         viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
%         viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
%         viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;
        
        viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
        viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
        viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
        viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;
        
        Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
        Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
        Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
        Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;
       
    
end
    coefmult=1;
    

% % %             for c5=2:levelsV
% % %             [Em1,En1,Ep1]=size(Edens{c5}.x2h);
% % %             viscmat{c5}=coefmult*restricthto2h3DV(viscmat{c5-1});
% % %             viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
% % %             viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
% % %             viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;
% % % %             
% % %             end


%make the coeficients for the relaxation steps in advance to save time in
%relaxation sweeps
% for c13=1:levelsP
% vcoef{c13}=[];
% pcoef{c13}=[];
% end


crossbarrier=zeros(numtimesteps,1);
for c13=1:levelsV
%     clear coef1 tEm tEn tEp

    
    tEm=Edens{c13}.Em;
    tEn=Edens{c13}.En;
    tEp=Edens{c13}.Ep;
    tEmm=tEm-1;
    tEmmm=tEm-2;
    tEnm=tEn-1;
    tEnmm=tEn-2;
    tEpm=tEp-1;
    tEpmm=tEp-2;
    

    vcoef{c13}.coefx=Re*st*Edens{c13}.Edensin(2:tEmm,2:tEnm,2:tEpm)/dt+(2*viscmatmid{c13}.lr(2:tEmm,2:tEnm,2:tEpm)+...
        2*viscmatmid{c13}.lr(2:tEmm,1:tEnmm,2:tEpm)+viscmatmid{c13}.ud(2:tEmm,2:tEnm,2:tEpm)+viscmatmid{c13}.ud(1:tEmmm,2:tEnm,2:tEpm)+...
        viscmatmid{c13}.fb(2:tEmm,2:tEnm,2:tEpm)+viscmatmid{c13}.fb(2:tEmm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    vcoef{c13}.coefxboundF=Re*st*Edens{c13}.Edensin(2:tEmm,2:tEnm,tEp)/dt+(2*viscmatmid{c13}.lr(2:tEmm,2:tEnm,tEp)+...
        2*viscmatmid{c13}.lr(2:tEmm,1:tEnmm,tEp)+viscmatmid{c13}.ud(2:tEmm,2:tEnm,tEp)+viscmatmid{c13}.ud(1:tEmmm,2:tEnm,tEp)+...
        viscmatmid{c13}.fb(2:tEmm,2:tEnm,1)+viscmatmid{c13}.fb(2:tEmm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    vcoef{c13}.coefxboundR=Re*st*Edens{c13}.Edensin(2:tEmm,tEn,2:tEpm)/dt+(2*viscmatmid{c13}.lr(2:tEmm,1,2:tEpm)+...
        2*viscmatmid{c13}.lr(2:tEmm,tEnm,2:tEpm)+viscmatmid{c13}.ud(2:tEmm,1,2:tEpm)+viscmatmid{c13}.ud(1:tEmmm,1,2:tEpm)+...
        viscmatmid{c13}.fb(2:tEmm,1,2:tEpm)+viscmatmid{c13}.fb(2:tEmm,1,1:tEpmm))/Edens{c13}.dx^2;
    vcoef{c13}.coefxboundFR=Re*st*Edens{c13}.Edensin(2:tEmm,tEn,tEp)/dt+(2*viscmatmid{c13}.lr(2:tEmm,1,tEp)+...
        2*viscmatmid{c13}.lr(2:tEmm,tEnm,tEp)+viscmatmid{c13}.ud(2:tEmm,1,tEp)+viscmatmid{c13}.ud(1:tEmmm,1,tEp)+...
        viscmatmid{c13}.fb(2:tEmm,1,1)+viscmatmid{c13}.fb(2:tEmm,1,tEpm))/Edens{c13}.dx^2;
    
    
    vcoef{c13}.coefy=Re*st*Edens{c13}.Edensin(2:tEmm,2:tEnm,2:tEpm)/dt+(viscmatmid{c13}.lr(2:tEmm,2:tEnm,2:tEpm)+...
        viscmatmid{c13}.lr(2:tEmm,1:tEnmm,2:tEpm)+2*viscmatmid{c13}.ud(2:tEmm,2:tEnm,2:tEpm)+2*viscmatmid{c13}.ud(1:tEmmm,2:tEnm,2:tEpm)+...
        viscmatmid{c13}.fb(2:tEmm,2:tEnm,2:tEpm)+viscmatmid{c13}.fb(2:tEmm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    vcoef{c13}.coefyboundF=Re*st*Edens{c13}.Edensin(2:tEmm,2:tEnm,tEp)/dt+(viscmatmid{c13}.lr(2:tEmm,2:tEnm,tEp)+...
        viscmatmid{c13}.lr(2:tEmm,1:tEnmm,tEp)+2*viscmatmid{c13}.ud(2:tEmm,2:tEnm,tEp)+2*viscmatmid{c13}.ud(1:tEmmm,2:tEnm,tEp)+...
        viscmatmid{c13}.fb(2:tEmm,2:tEnm,1)+viscmatmid{c13}.fb(2:tEmm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    vcoef{c13}.coefyboundR=Re*st*Edens{c13}.Edensin(2:tEmm,tEn,2:tEpm)/dt+(viscmatmid{c13}.lr(2:tEmm,1,2:tEpm)+...
        viscmatmid{c13}.lr(2:tEmm,tEnm,2:tEpm)+2*viscmatmid{c13}.ud(2:tEmm,1,2:tEpm)+2*viscmatmid{c13}.ud(1:tEmmm,1,2:tEpm)+...
        viscmatmid{c13}.fb(2:tEmm,1,2:tEpm)+viscmatmid{c13}.fb(2:tEmm,1,1:tEpmm))/Edens{c13}.dx^2;
    vcoef{c13}.coefyboundFR=Re*st*Edens{c13}.Edensin(2:tEmm,tEn,tEp)/dt+(viscmatmid{c13}.lr(2:tEmm,1,tEp)+...
        viscmatmid{c13}.lr(2:tEmm,tEnm,tEp)+2*viscmatmid{c13}.ud(2:tEmm,1,tEp)+2*viscmatmid{c13}.ud(1:tEmmm,1,tEp)+...
        viscmatmid{c13}.fb(2:tEmm,1,1)+viscmatmid{c13}.fb(2:tEmm,1,tEpm))/Edens{c13}.dx^2;
    
    
    vcoef{c13}.coefz=Re*st*Edens{c13}.Edensin(2:tEmm,2:tEnm,2:tEpm)/dt+(viscmatmid{c13}.lr(2:tEmm,2:tEnm,2:tEpm)+...
        viscmatmid{c13}.lr(2:tEmm,1:tEnmm,2:tEpm)+viscmatmid{c13}.ud(2:tEmm,2:tEnm,2:tEpm)+viscmatmid{c13}.ud(1:tEmmm,2:tEnm,2:tEpm)+...
        2*viscmatmid{c13}.fb(2:tEmm,2:tEnm,2:tEpm)+2*viscmatmid{c13}.fb(2:tEmm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    vcoef{c13}.coefzboundF=Re*st*Edens{c13}.Edensin(2:tEmm,2:tEnm,tEp)/dt+(viscmatmid{c13}.lr(2:tEmm,2:tEnm,tEp)+...
        viscmatmid{c13}.lr(2:tEmm,1:tEnmm,tEp)+viscmatmid{c13}.ud(2:tEmm,2:tEnm,tEp)+viscmatmid{c13}.ud(1:tEmmm,2:tEnm,tEp)+...
        2*viscmatmid{c13}.fb(2:tEmm,2:tEnm,1)+2*viscmatmid{c13}.fb(2:tEmm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    vcoef{c13}.coefzboundR=Re*st*Edens{c13}.Edensin(2:tEmm,tEn,2:tEpm)/dt+(viscmatmid{c13}.lr(2:tEmm,1,2:tEpm)+...
        viscmatmid{c13}.lr(2:tEmm,tEnm,2:tEpm)+viscmatmid{c13}.ud(2:tEmm,1,2:tEpm)+viscmatmid{c13}.ud(1:tEmmm,1,2:tEpm)+...
        2*viscmatmid{c13}.fb(2:tEmm,1,2:tEpm)+2*viscmatmid{c13}.fb(2:tEmm,1,1:tEpmm))/Edens{c13}.dx^2;
    vcoef{c13}.coefzboundFR=Re*st*Edens{c13}.Edensin(2:tEmm,tEn,tEp)/dt+(viscmatmid{c13}.lr(2:tEmm,1,tEp)+...
        viscmatmid{c13}.lr(2:tEmm,tEnm,tEp)+viscmatmid{c13}.ud(2:tEmm,1,tEp)+viscmatmid{c13}.ud(1:tEmmm,1,tEp)+...
        2*viscmatmid{c13}.fb(2:tEmm,1,1)+2*viscmatmid{c13}.fb(2:tEmm,1,tEpm))/Edens{c13}.dx^2;
    
    
    
    vcoef{c13}.rescoefpp1=viscmatmid{c13}.lr(2:tEmm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp1boundF=viscmatmid{c13}.lr(2:tEmm,2:tEnm,tEp)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp1boundR=viscmatmid{c13}.lr(2:tEmm,1,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp1boundFR=viscmatmid{c13}.lr(2:tEmm,1,tEp)/Edens{c13}.dx^2;
    
    vcoef{c13}.rescoefpm1=viscmatmid{c13}.lr(2:tEmm,1:tEnmm,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm1boundF=viscmatmid{c13}.lr(2:tEmm,1:tEnmm,tEp)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm1boundR=viscmatmid{c13}.lr(2:tEmm,tEnm,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm1boundFR=viscmatmid{c13}.lr(2:tEmm,tEnm,tEp)/Edens{c13}.dx^2;
    
    vcoef{c13}.rescoefpp2=viscmatmid{c13}.ud(2:tEmm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp2boundF=viscmatmid{c13}.ud(2:tEmm,2:tEnm,tEp)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp2boundR=viscmatmid{c13}.ud(2:tEmm,tEn,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp2boundFR=viscmatmid{c13}.ud(2:tEmm,tEn,tEp)/Edens{c13}.dx^2;
    
    vcoef{c13}.rescoefpm2=viscmatmid{c13}.ud(1:tEmmm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm2boundF=viscmatmid{c13}.ud(1:tEmmm,2:tEnm,tEp)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm2boundR=viscmatmid{c13}.ud(1:tEmmm,tEn,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm2boundFR=viscmatmid{c13}.ud(1:tEmmm,tEn,tEp)/Edens{c13}.dx^2;
    
    
    vcoef{c13}.rescoefpp3=viscmatmid{c13}.fb(2:tEmm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp3boundF=viscmatmid{c13}.fb(2:tEmm,2:tEnm,1)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp3boundR=viscmatmid{c13}.fb(2:tEmm,tEn,2:tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpp3boundFR=viscmatmid{c13}.fb(2:tEmm,tEn,1)/Edens{c13}.dx^2;
    
    vcoef{c13}.rescoefpm3=viscmatmid{c13}.fb(2:tEmm,2:tEnm,1:tEpmm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm3boundF=viscmatmid{c13}.fb(2:tEmm,2:tEnm,tEpm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm3boundR=viscmatmid{c13}.fb(2:tEmm,tEn,1:tEpmm)/Edens{c13}.dx^2;
    vcoef{c13}.rescoefpm3boundFR=viscmatmid{c13}.fb(2:tEmm,tEn,tEpm)/Edens{c13}.dx^2;
    
    
    vcoef{c13}.viscpp1=viscmat{c13}(2:tEmm,3:tEn,2:tEp)/4/Edens{c13}.dx^2;
    vcoef{c13}.viscpp1(:,tEnm,:)=viscmat{c13}(2:tEmm,2,2:tEp)/4/Edens{c13}.dx^2;
%     vcoef{c13}.viscpp1(:,tEnm,tEpm)=viscmat{c13}(2:tEmm,2,tEp)/4/Edens{c13}.dx^2;
    
    vcoef{c13}.viscpm1=viscmat{c13}(2:tEmm,1:tEnmm,2:tEp)/4/Edens{c13}.dx^2;
    vcoef{c13}.viscpm1(:,tEnm,:)=viscmat{c13}(2:tEmm,tEnm,2:tEp)/4/Edens{c13}.dx^2;
%     vcoef{c13}.viscpm1(:,tEnm,tEpm)=viscmat{c13}(2:tEmm,tEnm,tEpm)/4/Edens{c13}.dx^2;
    
    vcoef{c13}.viscpp2=viscmat{c13}(3:tEm,2:tEn,2:tEp)/4/Edens{c13}.dx^2;
    vcoef{c13}.viscpm2=viscmat{c13}(1:tEmmm,2:tEn,2:tEp)/4/Edens{c13}.dx^2;
    
    vcoef{c13}.viscpp3=viscmat{c13}(2:tEmm,2:tEn,3:tEp)/4/Edens{c13}.dx^2;
    vcoef{c13}.viscpp3(:,:,tEpm)=viscmat{c13}(2:tEmm,2:tEn,2)/4/Edens{c13}.dx^2; %add extra column for periodicity
%     vcoef{c13}.viscpp3(:,tEnm,tEpm)=viscmat{c13}(2:tEmm,2,2)/4/Edens{c13}.dx^2;
    
    vcoef{c13}.viscpm3=viscmat{c13}(2:tEmm,2:tEn,1:tEpmm)/4/Edens{c13}.dx^2;    
    vcoef{c13}.viscpm3(:,:,tEpm)=viscmat{c13}(2:tEmm,2:tEn,tEpm)/4/Edens{c13}.dx^2; % add extra column for periodicity
%     vcoef{c13}.viscpm3(:,tEnm,tEpm)=viscmat{c13}(2:tEmm,tEn,tEpm)/4/Edens{c13}.dx^2;        
       
    pcoef{c13}.coef=(1./Edens{c13}.Edensmidlr(2:tEmm,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidlr(2:tEmm,1:tEnmm,2:tEpm)+...
        1./Edens{c13}.Edensmidud(2:tEmm,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidud(1:tEmmm,2:tEnm,2:tEpm)+...
        1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    
    %Neumann boundary not at periodic ends
    pcoef{c13}.coefboundB=(1./Edens{c13}.Edensmidlr(1,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidlr(1,1:tEnmm,2:tEpm)+...
        2./Edens{c13}.Edensmidud(1,2:tEnm,2:tEpm)+...
        1./Edens{c13}.Edensmidfb(1,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidfb(1,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT=(1./Edens{c13}.Edensmidlr(tEm,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidlr(tEm,1:tEnmm,2:tEpm)+...
        2./Edens{c13}.Edensmidud(tEmm,2:tEnm,2:tEpm)+...
        1./Edens{c13}.Edensmidfb(tEm,2:tEnm,2:tEpm)+1./Edens{c13}.Edensmidfb(tEm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    
    %formula for non-edge points on periodic ends
    pcoef{c13}.coef(1:tEmmm,1:tEnmm,tEpm)=(1./Edens{c13}.Edensmidlr(2:tEmm,2:tEnm,1)+1./Edens{c13}.Edensmidlr(2:tEmm,1:tEnmm,1)+... %adjust indices to get right size
        1./Edens{c13}.Edensmidud(2:tEmm,2:tEnm,1)+1./Edens{c13}.Edensmidud(1:tEmmm,2:tEnm,1)+...
        1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,1)+1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coef(1:tEmmm,tEnm,1:tEpmm)=(1./Edens{c13}.Edensmidlr(2:tEmm,1,2:tEpm)+1./Edens{c13}.Edensmidlr(2:tEmm,tEnm,2:tEpm)+...
        1./Edens{c13}.Edensmidud(2:tEmm,1,2:tEpm)+1./Edens{c13}.Edensmidud(1:tEmmm,1,2:tEpm)+...
        1./Edens{c13}.Edensmidfb(2:tEmm,1,2:tEpm)+1./Edens{c13}.Edensmidfb(2:tEmm,1,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coef(1:tEmmm,tEnm,tEpm)=(1./Edens{c13}.Edensmidlr(2:tEmm,1,1)+1./Edens{c13}.Edensmidlr(2:tEmm,tEnm,1)+...
        1./Edens{c13}.Edensmidud(2:tEmm,1,1)+1./Edens{c13}.Edensmidud(1:tEmmm,1,1)+...
        1./Edens{c13}.Edensmidfb(2:tEmm,1,1)+1./Edens{c13}.Edensmidfb(2:tEmm,1,tEpm))/Edens{c13}.dx^2;
    
    %edge-point on periodic boundary
    pcoef{c13}.coefboundB(1,1:tEnmm,tEpm)=(1./Edens{c13}.Edensmidlr(1,2:tEnm,tEp)+1./Edens{c13}.Edensmidlr(1,1:tEnmm,tEp)+...
        2./Edens{c13}.Edensmidud(1,2:tEnm,tEp)+...
        1./Edens{c13}.Edensmidfb(1,2:tEnm,1)+1./Edens{c13}.Edensmidfb(1,2:tEnm,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundB(1,tEnm,1:tEpmm)=(1./Edens{c13}.Edensmidlr(1,1,2:tEpm)+1./Edens{c13}.Edensmidlr(1,tEnm,2:tEpm)+...
        2./Edens{c13}.Edensmidud(1,tEn,2:tEpm)+...
        1./Edens{c13}.Edensmidfb(1,tEn,2:tEpm)+1./Edens{c13}.Edensmidfb(1,tEn,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT(1,1:tEnmm,tEpm)=(1./Edens{c13}.Edensmidlr(tEm,2:tEnm,tEp)+1./Edens{c13}.Edensmidlr(tEm,1:tEnmm,tEp)+...
        2./Edens{c13}.Edensmidud(tEmm,2:tEnm,tEp)+...
        1./Edens{c13}.Edensmidfb(tEm,2:tEnm,1)+1./Edens{c13}.Edensmidfb(tEm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT(1,tEnm,1:tEpmm)=(1./Edens{c13}.Edensmidlr(tEm,1,2:tEpm)+1./Edens{c13}.Edensmidlr(tEm,tEnm,2:tEpm)+...
        2./Edens{c13}.Edensmidud(tEmm,tEn,2:tEpm)+...
        1./Edens{c13}.Edensmidfb(tEm,tEn,2:tEpm)+1./Edens{c13}.Edensmidfb(tEm,tEn,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundB(1,tEnm,tEpm)=(1./Edens{c13}.Edensmidlr(1,1,tEp)+1./Edens{c13}.Edensmidlr(1,tEnm,tEp)+...
        2./Edens{c13}.Edensmidud(1,tEn,tEp)+...
        1./Edens{c13}.Edensmidfb(1,tEn,1)+1./Edens{c13}.Edensmidfb(1,tEn,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT(1,tEnm,tEpm)=(1./Edens{c13}.Edensmidlr(tEm,1,tEp)+1./Edens{c13}.Edensmidlr(tEm,tEnm,tEp)+...
        2./Edens{c13}.Edensmidud(tEmm,tEn,tEp)+...
        1./Edens{c13}.Edensmidfb(tEm,tEn,1)+1./Edens{c13}.Edensmidfb(tEm,tEn,tEpm))/Edens{c13}.dx^2;

       
    pcoef{c13}.rescoefpp1=1./Edens{c13}.Edensmidlr(2:tEmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp1(1:tEmmm,tEnm,1:tEpm)=1./Edens{c13}.Edensmidlr(2:tEmm,1,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp1boundB=1./Edens{c13}.Edensmidlr(1,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp1boundB(1,tEnm,1:tEpm)=1./Edens{c13}.Edensmidlr(1,1,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp1boundT=1./Edens{c13}.Edensmidlr(tEm,2:tEnm,2:tEp)/Edens{c13}.dx^2; %
    pcoef{c13}.rescoefpp1boundT(1,tEnm,1:tEpm)=1./Edens{c13}.Edensmidlr(tEm,1,2:tEp)/Edens{c13}.dx^2;
 
    pcoef{c13}.rescoefpm1=1./Edens{c13}.Edensmidlr(2:tEmm,1:tEnmm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1(1:tEmmm,tEnm,1:tEpm)=1./Edens{c13}.Edensmidlr(2:tEmm,tEnm,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpm1boundB=1./Edens{c13}.Edensmidlr(1,1:tEnmm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1boundB(1,tEnm,1:tEpm)=1./Edens{c13}.Edensmidlr(1,tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1boundT=1./Edens{c13}.Edensmidlr(tEm,1:tEnmm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1boundT(1,tEnm,1:tEpm)=1./Edens{c13}.Edensmidlr(tEm,tEnm,2:tEp)/Edens{c13}.dx^2;


    pcoef{c13}.rescoefpp2=1./Edens{c13}.Edensmidud(2:tEmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp2(1:tEmmm,tEnm,1:tEpm)=1./Edens{c13}.Edensmidud(2:tEmm,tEn,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp2boundB=2./Edens{c13}.Edensmidud(1,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp2boundB(1,tEnm,1:tEpm)=2./Edens{c13}.Edensmidud(1,tEn,2:tEp)/Edens{c13}.dx^2;

    pcoef{c13}.rescoefpm2=1./Edens{c13}.Edensmidud(1:tEmmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm2(1:tEmmm,tEnm,1:tEpm)=1./Edens{c13}.Edensmidud(1:tEmmm,tEn,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpm2boundT=2./Edens{c13}.Edensmidud(tEmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm2boundT(1,tEnm,1:tEpm)=2./Edens{c13}.Edensmidud(tEmm,tEn,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp3=1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3(1:tEmmm,1:tEnmm,tEpm)=1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,1)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3(1:tEmmm,tEnm,1:tEpmm)=1./Edens{c13}.Edensmidfb(2:tEmm,tEn,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3(1:tEmmm,tEnm,tEpm)=1./Edens{c13}.Edensmidfb(2:tEmm,tEn,1)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp3boundB=1./Edens{c13}.Edensmidfb(1,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundB(1,1:tEnmm,tEpm)=1./Edens{c13}.Edensmidfb(1,2:tEnm,1)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundB(1,tEnm,1:tEpmm)=1./Edens{c13}.Edensmidfb(1,tEn,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundB(1,tEnm,tEpm)=1./Edens{c13}.Edensmidfb(1,tEn,1)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp3boundT=1./Edens{c13}.Edensmidfb(tEm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundT(1,1:tEnmm,tEpm)=1./Edens{c13}.Edensmidfb(tEm,2:tEnm,1)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundT(1,tEnm,1:tEpmm)=1./Edens{c13}.Edensmidfb(tEm,tEn,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundT(1,tEnm,tEpm)=1./Edens{c13}.Edensmidfb(tEm,tEn,1)/Edens{c13}.dx^2;

    pcoef{c13}.rescoefpm3=1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,1:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3(1:tEmmm,1:tEnmm,tEpm)=1./Edens{c13}.Edensmidfb(2:tEmm,2:tEnm,tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3(1:tEmmm,tEnm,1:tEpmm)=1./Edens{c13}.Edensmidfb(2:tEmm,tEn,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3(1:tEmmm,tEnm,tEpm)=1./Edens{c13}.Edensmidfb(2:tEmm,tEn,tEpm)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpm3boundB=1./Edens{c13}.Edensmidfb(1,2:tEnm,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundB(1,:,tEpm)=1./Edens{c13}.Edensmidfb(1,2:tEnm,tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundB(1,tEnm,1:tEpmm)=1./Edens{c13}.Edensmidfb(1,tEn,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundB(1,tEnm,tEpm)=1./Edens{c13}.Edensmidfb(1,tEn,tEpm)/Edens{c13}.dx^2;

    pcoef{c13}.rescoefpm3boundT=1./Edens{c13}.Edensmidfb(tEm,2:tEnm,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundT(1,:,tEpm)=1./Edens{c13}.Edensmidfb(tEm,2:tEnm,tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundT(1,tEnm,1:tEpmm)=1./Edens{c13}.Edensmidfb(tEm,tEn,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundT(1,tEnm,tEpm)=1./Edens{c13}.Edensmidfb(tEm,tEn,tEpm)/Edens{c13}.dx^2;
    
    
    pcoef{c13}.coefpp1=pcoef{c13}.rescoefpp1./pcoef{c13}.coef;
    pcoef{c13}.coefpp1boundT=pcoef{c13}.rescoefpp1boundT./pcoef{c13}.coefboundT;
    pcoef{c13}.coefpp1boundB=pcoef{c13}.rescoefpp1boundB./pcoef{c13}.coefboundB;
    
    pcoef{c13}.coefpm1=pcoef{c13}.rescoefpm1./pcoef{c13}.coef;
    pcoef{c13}.coefpm1boundB=pcoef{c13}.rescoefpm1boundB./pcoef{c13}.coefboundB;
    pcoef{c13}.coefpm1boundT=pcoef{c13}.rescoefpm1boundT./pcoef{c13}.coefboundT;
    
    
    pcoef{c13}.coefpp2=pcoef{c13}.rescoefpp2./pcoef{c13}.coef;
    pcoef{c13}.coefpp2boundB=pcoef{c13}.rescoefpp2boundB./pcoef{c13}.coefboundB;
    
    pcoef{c13}.coefpm2=pcoef{c13}.rescoefpm2./pcoef{c13}.coef;
    pcoef{c13}.coefpm2boundT=pcoef{c13}.rescoefpm2boundT./pcoef{c13}.coefboundT;
        
    pcoef{c13}.coefpp3=pcoef{c13}.rescoefpp3./pcoef{c13}.coef;
    pcoef{c13}.coefpp3boundB=pcoef{c13}.rescoefpp3boundB./pcoef{c13}.coefboundB;
    pcoef{c13}.coefpp3boundT=pcoef{c13}.rescoefpp3boundT./pcoef{c13}.coefboundT;
    
    pcoef{c13}.coefpm3=pcoef{c13}.rescoefpm3./pcoef{c13}.coef;
    pcoef{c13}.coefpm3boundB=pcoef{c13}.rescoefpm3boundB./pcoef{c13}.coefboundB;
    pcoef{c13}.coefpm3boundT=pcoef{c13}.rescoefpm3boundT./pcoef{c13}.coefboundT;
   

end



  
% max(max(Edens))
% min(min(Edens))  X(z,x,y)
stuckb=find(X(:,2)<10^-15);
stuckt=find(X(:,2)>ylength-10^-15);
Xpast=zeros(numOfPoints,1);
%% Go through all of the timesteps
uxsave=ux;
uysave=uy;
uzsave=uz;
psave=p;
vsave=viscmat{1};
densave=Edens{1}.Edensin;

issue=1;

%initialize tracers
S=zeros(108,3);
for i=1:6
    for j=1:6
        S(6*(i-1)+j,:)=[i*xlength/7, ylength, j*zlength/7];
    end
end

S(37:72,:)=S(1:36,:);
S(73:108,:)=S(1:36,:);
S(37:72,2)=S(37:72,2)-h;
S(73:108,2)=S(73:108,2)-2*h;

SSave=S;
S0=S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for c12=1:levels
%     Ind{c12}.boundary=find(Ind{1}.boundary);
%     Ind{c12}.boundaryL=find(Ind{1}.boundaryL);
%     Ind{c12}.boundaryR=find(Ind{1}.boundaryR);
%     Ind{c12}.interior=find(Ind{1}.interior);
%     Ind{c12}.interiorL=find(Ind{1}.interiorL);
%     Ind{c12}.interiorR=find(Ind{1}.interiorR);
%     Ind{c12}.interiorU=find(Ind{1}.interiorU);
%     Ind{c12}.interiorD=find(Ind{1}.interiorD);
%     Ind{c12}.interiorF=find(Ind{1}.interiorF);
%     Ind{c12}.interiorB=find(Ind{1}.interiorB);
%     Ind{c12}.exterior=find(Ind{1}.exterior);
% end
eShear=zeros(numtimesteps,1); vShear=eShear; Shear_Force=eShear;
perr=zeros(numtimesteps,1); verr=perr;


Enm=En-1;
Enmm=En-2;
Emm=Em-1;
Emmm=Em-2;
uhalfx=ux;
uhalfy=uy;
uhalfz=uz;
Emm=Em-1;
Enm=En-1;
Emmm=Em-2;
Enmm=En-2;
Epm=Ep-1;
Epmm=Ep-2;
% matlabpool(levels) 

% fcd0mean=fc/d0mean^3;
% matlabpool(3)

t_star=floor(pi/w/dt);
e0=e0/charLength;

Strain=e0*sin([0:numtimesteps]*dt).*(1./(1+exp(-2*[0:numtimesteps]*dt))-1/2).^2;

fStrain=zeros(numtimesteps,1);
bStrain=fStrain;
tStrain=fStrain;
bStraintop=bStrain;
fStraintop=fStrain;
