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
p0=1;

%%%%%%%%%%%%nondimensional constants
st=timefreq*charLength/v0;
eu=p0/rho0/v0^2;
Re=rho0*charLength*v0/visc0;
% fmax=1;
% f0=175000;%3 is because 3 times rest length is the breaking criteria
fc=charLength/rho0/v0^2;
B=b;
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
uxs=ux; uys=uy; uzs=uz;

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

b=b/e0./d0;
b(~A)=0;
t=0:dt:tend;
[rowind,colind]=find(A);
matind=find(A);

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

[Edens,viscmat, viscmatmid]=dens_visc(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,ylength, xlength,visc,addlvisc,levelsV,Em,En,Ep);
coefmult=1;
S{1}=x;
S{2}=y;
S{3}=z;

crossbarrier=zeros(numtimesteps,1);
[vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
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
Ssave=cell(3,1);
Ssave{1}=x;
Ssave{2}=y;
Ssave{3}=z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

Break=0;
c3=1;
uxm=ux;
uym=uy;
uzm=uz;
ux1=ux;
ux2=ux;
ux3=ux;
uy1=uy;
uy2=uy;
uy3=uy;
uz1=uz;
uz2=uz;
uz3=uz;


dx=h;
   
    
uz(1,:,:)=0;
uz(Em,:,:)=(exp(2*c3*dt)-1)*((exp(4*c3*dt)-1)*cos(c3*dt)+exp(2*c3*dt)*8*sin(c3*dt))/(4*(1+exp(2*c3*dt))^3);
uzs(1,:,:)=uz(1,:,:);
uzs(Em,:,:)=uz(Em,:,:);
[uhalfx,uhalfy,uhalfz,verr(c3)]=multigVEL3D_CGper2(dt,h,ux,uy,uz,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx*fc,Efy*fc,Efz*fc,...
    Em,En,Ep,levelsV,vcoef,st,eust,Re);

%Now solve for new pressure
[p,perr(c3)]=multigPRESSUREprod3Dper2(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levelsP,pcoef,eust);

%Solve for the new velocity profile u

edenscoef=eust*dt/(2*h)./Edens{1}.Edensin(2:Emm,2:En,2:Ep);

ux(2:Emm,2:Enm,2:Epm)=(p(2:Emm,1:Enmm,2:Epm)-...
    p(2:Emm,3:En,2:Epm)).*edenscoef(:,1:Enmm,1:Epmm)+uhalfx(2:Emm,2:Enm,2:Epm);

ux(2:Emm,2:Enm,Ep)=(p(2:Emm,1:Enmm,Ep)-...
    p(2:Emm,3:En,Ep)).*edenscoef(:,1:Enmm,Epm)+uhalfx(2:Emm,2:Enm,Ep);
ux(2:Emm,2:Enm,1)=ux(2:Emm,2:Enm,Ep);

ux(2:Emm,En,2:Epm)=(p(2:Emm,Enm,2:Epm)-...
    p(2:Emm,2,2:Epm)).*edenscoef(:,Enm,1:Epmm)+uhalfx(2:Emm,En,2:Epm);
ux(2:Emm,1,2:Epm)=ux(2:Emm,En,2:Epm);

ux(2:Emm,En,Ep)=(p(2:Emm,Enm,Ep)-p(2:Emm,2,Ep)).*edenscoef(:,Enm,Epm)+uhalfx(2:Emm,En,Ep);
ux(2:Emm,En,1)=ux(2:Emm,En,Ep);
ux(2:Emm,1,Ep)=ux(2:Emm,En,Ep);
ux(2:Emm,1,1)=ux(2:Emm,En,Ep);


uy(2:Emm,2:En,2:Epm)=(p(1:Emmm,2:En,2:Epm)-...
    p(3:Em,2:En,2:Epm)).*edenscoef(:,:,1:Epmm)+uhalfy(2:Emm,2:En,2:Epm);
uy(2:Emm,2:En,Ep)=(p(1:Emmm,2:En,Ep)-...
    p(3:Em,2:En,Ep)).*edenscoef(:,:,Epm)+uhalfy(2:Emm,2:En,Ep);
uy(2:Emm,2:En,1)=uy(2:Emm,2:En,Ep);
uy(2:Emm,1,:)=uy(2:Emm,En,:);

uz(2:Emm,2:En,2:Epm)=(p(2:Emm,2:En,1:Epmm)-...
    p(2:Emm,2:En,3:Ep)).*edenscoef(:,:,1:Epmm)+uhalfz(2:Emm,2:En,2:Epm);
uz(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Epm)-...
    p(2:Emm,2:En,2)).*edenscoef(:,:,Epm)+uhalfz(2:Emm,2:En,Ep);
uz(2:Emm,2:En,1)=uz(2:Emm,2:En,Ep);
uz(2:Emm,1,:)=uz(2:Emm,En,:);
 if addlvisc>0 || addldens>0
    viscmat{1}=transferLtoEvisc3Dper_e2(h,x,y,z,X,visc,addlvisc,numOfnonzero,1,zlength,xlength);
    viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
    viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
    viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;

     Edens{1}.Edensin=transferLtoEdens3Dper_e2(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);

    for c5=2:levelsV
        [Em1,En1,Ep1]=size(Edens{c5}.x2h);


        viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
        viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
        viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
        viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;

        Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
        Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
        Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
        Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;

    end
 end
 
 
 [vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);



for c3=1:numtimesteps   
    
    c3
     U1=transferEtoLvel3Dper_e2a(h,ux,uy,uz,x,y,z,X,zlength,xlength); %efficient E-L transfer function
    if isnan(mean(U(:,3)))==1 || isnan(mean(U(:,2)))==1 ||isnan(mean(U(:,1)))==1 
       'error, U not computed'
       break
    end
    
    U1(stuckb,:)=0;
    U1(stuckt,3)=uz(Em,2,3); U1(stuckt,1:2)=0;
    %calculate new position from the new velocity using Euler's Method
    Xback=zeros(numOfPoints,1);
    X1=U1*dt/2/st+X;
        
    Xpast(X1(:,3)>zlength)=Xpast(X1(:,3)>zlength)+zlength;
    Xback(X1(:,3)<0)=Xback(X1(:,3)<0)-zlength;
    X1(:,3)=mod(X1(:,3),zlength);
    X1(:,1)=mod(X1(:,1),xlength);
    X1(stuckb,2)=0;
    X1(stuckt,2)=ylength;
    X1((X1(:,2)>ylength),2)=ylength;
    X1((X1(:,2)<0),2)=0;
    
    X1(lower,2)=0;
    X1(upper,2)=ylength;
    
    
    Zdistr=zeros(size(Zdist));
    Zdistl=Zdistr;
    Xdistr=Zdistr; 
    Xdistl=Zdistl;
    %calculate forces at the new positions
    
%     Xcalc=[X(:,1),X(:,2),X(:,3)+zlength*Xpast];
    Xdist(matind)=X1(colind,1)-X1(rowind,1);
    Xdistr(matind)=X1(colind,1)-X1(rowind,1)+xlength;
    Xdistl(matind)=X1(colind,1)-X1(rowind,1)-xlength;
    Ydist(matind)=X1(colind,2)-X1(rowind,2);
    Zdist(matind)=X1(colind,3)-X1(rowind,3);
    Zdistr(matind)=X1(colind,3)-X1(rowind,3)+zlength;
    Zdistl(matind)=X1(colind,3)-X1(rowind,3)-zlength;
    
    zind=find(abs(Zdist)>abs(Zdistr));
    Zdist(zind)=Zdistr(zind);
    zind=find(abs(Zdist)>abs(Zdistl));
    Zdist(zind)=Zdistl(zind);
    
    xind=find(abs(Xdist)>abs(Xdistr));
    Xdist(xind)=Xdistr(xind);
    xind=find(abs(Xdist)>abs(Xdistl));
    Xdist(xind)=Xdistl(xind);
    
    DSq=Xdist.^2+Ydist.^2+Zdist.^2;
    D=sqrt(DSq);
    
    breakind=find(D>2*d0);
    if ~isempty(breakind)

        A(breakind)=0;
        [rowind,colind]=find(A);
        matind=find(A);
        Break=Break+1;
    end
    Dtemp=(D-d0).*K./D;
    Fx=Xdist.*Dtemp;
    Fy=Ydist.*Dtemp;
    Fz=Zdist.*Dtemp;
    
    Fx(~A)=0;
    Fy(~A)=0;
    Fz(~A)=0;
    Fx(isnan(Fx))=0;
    Fy(isnan(Fy))=0;
    Fz(isnan(Fz))=0;
    XrelVel(matind)=v0*(U(colind,1)-U(rowind,1));
    YrelVel(matind)=v0*(U(colind,2)-U(rowind,2));
    ZrelVel(matind)=v0*(U(colind,3)-U(rowind,3));    
    
    VrelDotXrel=b.*(XrelVel.*Xdist+YrelVel.*Ydist+ZrelVel.*Zdist)./DSq;
    XFb=VrelDotXrel.*Xdist;
    YFb=VrelDotXrel.*Ydist;
    ZFb=VrelDotXrel.*Zdist;

%     XFb=b*XrelVel;
%     YFb=b*YrelVel;
%     ZFb=b*ZrelVel;

    XFb(~A)=0;
    YFb(~A)=0;
    ZFb(~A)=0;
    XFb(isnan(XFb))=0;
    YFb(isnan(YFb))=0;
    ZFb(isnan(ZFb))=0;
    
    Fx(stuckt,:)=0;
    Fy(stuckt,:)=0;
    Fz(stuckt,:)=0;
    XFb(stuckb,:)=0;
    YFb(stuckb,:)=0;
    ZFb(stuckb,:)=0;

    Fx(stuckb,:)=0;
    Fy(stuckb,:)=0;
    Fz(stuckb,:)=0;
    XFb(stuckt,:)=0;
    YFb(stuckt,:)=0;
    ZFb(stuckt,:)=0;
    
    %Transfer the new forces and densities to the Eulerian points
    [Efx1,Efy1,Efz1]=transferLtoE3Dper_e2(h,x,y,z,X1,sum(Fx,2)+sum(XFb,2),sum(Fy,2)+sum(YFb,2),sum(Fz,2)+sum(ZFb,2),d0mean,zlength,xlength);
    
    
    if addlvisc>0 || addldens>0
        viscmat{1}=transferLtoEvisc3Dper_e2(h,x,y,z,X1,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
        
     Edens{1}.Edensin=transferLtoEdens3Dper_e2(h,x,y,z,X1,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);

        for c5=2:levelsV
            [Em1,En1,Ep1]=size(Edens{c5}.x2h);
            
            
            viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
            viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
            viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
            viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;

            Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
            Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
            Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
            Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;

        end
    end
        
        
  [vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt/2);


    uxm=ux;
    uym=uy;
    uzm=uz;
    
    dx=h;
    tstart=tic;
   
    
    uz(1,:,:)=0;
    uz(Em,:,:)=(exp(2*c3*dt)-1)*((exp(4*c3*dt)-1)*cos(c3*dt)+exp(2*c3*dt)*8*sin(c3*dt))/(4*(1+exp(2*c3*dt))^3);
    uz1(1,:,:)=uz(1,:,:);
    uz1(Em,:,:)=uz(Em,:,:);
    uz2(1,:,:)=uz(1,:,:);
    uz2(Em,:,:)=uz(Em,:,:);
    uz3(1,:,:)=uz(1,:,:);
    uz3(Em,:,:)=uz(Em,:,:);


[uhalfx,uhalfy,uhalfz,verr(c3)]=multigVEL3D_CGper2(dt/2,h,ux,uy,uz,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx1*fc,Efy1*fc,Efz1*fc,...
    Em,En,Ep,levelsV,vcoef,st,eust,Re);

    if size(p)~=[Em,En,Ep];
        p=zeros(Em,En,Ep);
        issue=1+issue;
    end
     if size(uhalfx)~=[Em,En,Ep];
        uxhalf=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfy)~=[Em,En,Ep];
        uhalfy=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfz)~=[Em,En,Ep];
        uhalfz=zeros(Em,En,Ep);
        issue=1+issue;
     end
        
    %Now solve for new pressure
    [p,perr(c3)]=multigPRESSUREprod3Dper2(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levelsP,pcoef,eust);
   
    %Solve for the new velocity profile u
    
%     parfor 


edenscoef=eust*dt/2/(2*h)./Edens{1}.Edensin(2:Emm,2:En,2:Ep);
  
    ux1(2:Emm,2:Enm,2:Epm)=(p(2:Emm,1:Enmm,2:Epm)-...
        p(2:Emm,3:En,2:Epm)).*edenscoef(:,1:Enmm,1:Epmm)+uhalfx(2:Emm,2:Enm,2:Epm);
    
    ux1(2:Emm,2:Enm,Ep)=(p(2:Emm,1:Enmm,Ep)-...
        p(2:Emm,3:En,Ep)).*edenscoef(:,1:Enmm,Epm)+uhalfx(2:Emm,2:Enm,Ep);
    ux1(2:Emm,2:Enm,1)=ux1(2:Emm,2:Enm,Ep);
   
    ux1(2:Emm,En,2:Epm)=(p(2:Emm,Enm,2:Epm)-...
        p(2:Emm,2,2:Epm)).*edenscoef(:,Enm,1:Epmm)+uhalfx(2:Emm,En,2:Epm);
    ux1(2:Emm,1,2:Epm)=ux1(2:Emm,En,2:Epm);
    
    ux1(2:Emm,En,Ep)=(p(2:Emm,Enm,Ep)-p(2:Emm,2,Ep)).*edenscoef(:,Enm,Epm)+uhalfx(2:Emm,En,Ep);
    ux1(2:Emm,En,1)=ux1(2:Emm,En,Ep);
    ux1(2:Emm,1,Ep)=ux1(2:Emm,En,Ep);
    ux1(2:Emm,1,1)=ux1(2:Emm,En,Ep);
    
    
    uy1(2:Emm,2:En,2:Epm)=(p(1:Emmm,2:En,2:Epm)-...
        p(3:Em,2:En,2:Epm)).*edenscoef(:,:,1:Epmm)+uhalfy(2:Emm,2:En,2:Epm);
    uy1(2:Emm,2:En,Ep)=(p(1:Emmm,2:En,Ep)-...
        p(3:Em,2:En,Ep)).*edenscoef(:,:,Epm)+uhalfy(2:Emm,2:En,Ep);
    uy1(2:Emm,2:En,1)=uy1(2:Emm,2:En,Ep);
    uy1(2:Emm,1,:)=uy1(2:Emm,En,:);
    
    uz1(2:Emm,2:En,2:Epm)=(p(2:Emm,2:En,1:Epmm)-...
        p(2:Emm,2:En,3:Ep)).*edenscoef(:,:,1:Epmm)+uhalfz(2:Emm,2:En,2:Epm);
    uz1(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Epm)-...
        p(2:Emm,2:En,2)).*edenscoef(:,:,Epm)+uhalfz(2:Emm,2:En,Ep);
    uz1(2:Emm,2:En,1)=uz1(2:Emm,2:En,Ep);
    uz1(2:Emm,1,:)=uz1(2:Emm,En,:);
    

    %Transfer the velocity to the Lagrangian points
    Um=U;
    U2=transferEtoLvel3Dper_e2a(h,ux1,uy1,uz1,x,y,z,X1,zlength,xlength); %efficient E-L transfer function
    if isnan(mean(U2(:,3)))==1 || isnan(mean(U2(:,2)))==1 ||isnan(mean(U2(:,1)))==1 
       'error, Us not computed'
       break   
    end
    
    U2(stuckb,:)=0;
    U2(stuckt,3)=uz1(Em,2,3); U2(stuckt,1:2)=0;
    %calculate new position from the new velocity using Euler's Method
    Xback=zeros(numOfPoints,1);
    X2=U2*dt/2/st+X;
        
    Xpast(X2(:,3)>zlength)=Xpast(X2(:,3)>zlength)+zlength;
    Xback(X2(:,3)<0)=Xback(X2(:,3)<0)-zlength;
    X2(:,3)=mod(X2(:,3),zlength);
    X2(:,1)=mod(X2(:,1),xlength);
    X2(stuckb,2)=0;
    X2(stuckt,2)=ylength;
    X2((X2(:,2)>ylength),2)=ylength;
    X2((X2(:,2)<0),2)=0;
    
    X2(lower,2)=0;
    X2(upper,2)=ylength;
    
    
    Zdistr=zeros(size(Zdist));
    Zdistl=Zdistr;
    Xdistr=Zdistr; 
    Xdistl=Zdistl;
    %calculate forces at the new positions
    
%     Xcalc=[X(:,1),X(:,2),X(:,3)+zlength*Xpast];
    Xdist(matind)=X2(colind,1)-X2(rowind,1);
    Xdistr(matind)=X2(colind,1)-X2(rowind,1)+xlength;
    Xdistl(matind)=X2(colind,1)-X2(rowind,1)-xlength;
    Ydist(matind)=X2(colind,2)-X2(rowind,2);
    Zdist(matind)=X2(colind,3)-X2(rowind,3);
    Zdistr(matind)=X2(colind,3)-X2(rowind,3)+zlength;
    Zdistl(matind)=X2(colind,3)-X2(rowind,3)-zlength;
    
    zind=find(abs(Zdist)>abs(Zdistr));
    Zdist(zind)=Zdistr(zind);
    zind=find(abs(Zdist)>abs(Zdistl));
    Zdist(zind)=Zdistl(zind);
    
    xind=find(abs(Xdist)>abs(Xdistr));
    Xdist(xind)=Xdistr(xind);
    xind=find(abs(Xdist)>abs(Xdistl));
    Xdist(xind)=Xdistl(xind);
    
    DSq=Xdist.^2+Ydist.^2+Zdist.^2;
    D=sqrt(DSq);
    
    breakind=find(D>2*d0);
    if ~isempty(breakind)

        A(breakind)=0;
        [rowind,colind]=find(A);
        matind=find(A);
        Break=Break+1;
    end
    Dtemp=(D-d0).*K./D;
    Fx=Xdist.*Dtemp;
    Fy=Ydist.*Dtemp;
    Fz=Zdist.*Dtemp;
    
    Fx(~A)=0;
    Fy(~A)=0;
    Fz(~A)=0;
    Fx(isnan(Fx))=0;
    Fy(isnan(Fy))=0;
    Fz(isnan(Fz))=0;
    XrelVel(matind)=v0*(U2(colind,1)-U2(rowind,1));
    YrelVel(matind)=v0*(U2(colind,2)-U2(rowind,2));
    ZrelVel(matind)=v0*(U2(colind,3)-U2(rowind,3));    
    
    VrelDotXrel=b.*(XrelVel.*Xdist+YrelVel.*Ydist+ZrelVel.*Zdist)./DSq;
    XFb=VrelDotXrel.*Xdist;
    YFb=VrelDotXrel.*Ydist;
    ZFb=VrelDotXrel.*Zdist;

%     XFb=b*XrelVel;
%     YFb=b*YrelVel;
%     ZFb=b*ZrelVel;

    XFb(~A)=0;
    YFb(~A)=0;
    ZFb(~A)=0;
    XFb(isnan(XFb))=0;
    YFb(isnan(YFb))=0;
    ZFb(isnan(ZFb))=0;
   
    Fx(stuckt,:)=0;
    Fy(stuckt,:)=0;
    Fz(stuckt,:)=0;
    XFb(stuckb,:)=0;
    YFb(stuckb,:)=0;
    ZFb(stuckb,:)=0;

    Fx(stuckb,:)=0;
    Fy(stuckb,:)=0;
    Fz(stuckb,:)=0;
    XFb(stuckt,:)=0;
    YFb(stuckt,:)=0;
    ZFb(stuckt,:)=0;
    
    %Transfer the new forces and densities to the Eulerian points
    [Efx2,Efy2,Efz2]=transferLtoE3Dper_e2(h,x,y,z,X2,sum(Fx,2)+sum(XFb,2),sum(Fy,2)+sum(YFb,2),sum(Fz,2)+sum(ZFb,2),d0mean,zlength,xlength);


    if addlvisc>0 || addldens>0
        viscmat{1}=transferLtoEvisc3Dper_e2(h,x,y,z,X2,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
        
     Edens{1}.Edensin=transferLtoEdens3Dper_e2(h,x,y,z,X2,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);

        for c5=2:levelsV
            [Em1,En1,Ep1]=size(Edens{c5}.x2h);
            
            
            viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
            viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
            viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
            viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;

            Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
            Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
            Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
            Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;
       


        end
    end
        
        
[vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
    

[uhalfx,uhalfy,uhalfz,verr(c3)]=multigVEL3D_CGper2(dt,h,ux1,uy1,uz1,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx2*fc,Efy2*fc,Efz2*fc,...
    Em,En,Ep,levelsV,vcoef,st,eust,Re);

    if size(p)~=[Em,En,Ep];
        p=zeros(Em,En,Ep);
        issue=1+issue;
    end
     if size(uhalfx)~=[Em,En,Ep];
        uxhalf=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfy)~=[Em,En,Ep];
        uhalfy=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfz)~=[Em,En,Ep];
        uhalfz=zeros(Em,En,Ep);
        issue=1+issue;
     end
        
    %Now solve for new pressure
    [p,perr(c3)]=multigPRESSUREprod3Dper2(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levelsP,pcoef,eust);
   
    %Solve for the new velocity profile u
    
%     parfor 


edenscoef=eust*dt/(2*h)./Edens{1}.Edensin(2:Emm,2:En,2:Ep);
 
    ux2(2:Emm,2:Enm,2:Epm)=(p(2:Emm,1:Enmm,2:Epm)-...
        p(2:Emm,3:En,2:Epm)).*edenscoef(:,1:Enmm,1:Epmm)+uhalfx(2:Emm,2:Enm,2:Epm);
    
    ux2(2:Emm,2:Enm,Ep)=(p(2:Emm,1:Enmm,Ep)-...
        p(2:Emm,3:En,Ep)).*edenscoef(:,1:Enmm,Epm)+uhalfx(2:Emm,2:Enm,Ep);
    ux2(2:Emm,2:Enm,1)=ux2(2:Emm,2:Enm,Ep);
   
    ux2(2:Emm,En,2:Epm)=(p(2:Emm,Enm,2:Epm)-...
        p(2:Emm,2,2:Epm)).*edenscoef(:,Enm,1:Epmm)+uhalfx(2:Emm,En,2:Epm);
    ux2(2:Emm,1,2:Epm)=ux2(2:Emm,En,2:Epm);
    
    ux2(2:Emm,En,Ep)=(p(2:Emm,Enm,Ep)-p(2:Emm,2,Ep)).*edenscoef(:,Enm,Epm)+uhalfx(2:Emm,En,Ep);
    ux2(2:Emm,En,1)=ux2(2:Emm,En,Ep);
    ux2(2:Emm,1,Ep)=ux2(2:Emm,En,Ep);
    ux2(2:Emm,1,1)=ux2(2:Emm,En,Ep);
    
    
    uy2(2:Emm,2:En,2:Epm)=(p(1:Emmm,2:En,2:Epm)-...
        p(3:Em,2:En,2:Epm)).*edenscoef(:,:,1:Epmm)+uhalfy(2:Emm,2:En,2:Epm);
    uy2(2:Emm,2:En,Ep)=(p(1:Emmm,2:En,Ep)-...
        p(3:Em,2:En,Ep)).*edenscoef(:,:,Epm)+uhalfy(2:Emm,2:En,Ep);
    uy2(2:Emm,2:En,1)=uy2(2:Emm,2:En,Ep);
    uy2(2:Emm,1,:)=uy2(2:Emm,En,:);
    
    uz2(2:Emm,2:En,2:Epm)=(p(2:Emm,2:En,1:Epmm)-...
        p(2:Emm,2:En,3:Ep)).*edenscoef(:,:,1:Epmm)+uhalfz(2:Emm,2:En,2:Epm);
    uz2(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Epm)-...
        p(2:Emm,2:En,2)).*edenscoef(:,:,Epm)+uhalfz(2:Emm,2:En,Ep);
    uz2(2:Emm,2:En,1)=uz2(2:Emm,2:En,Ep);
    uz2(2:Emm,1,:)=uz2(2:Emm,En,:);
    
  







  U3=transferEtoLvel3Dper_e2a(h,ux2,uy2,uz2,x,y,z,X2,zlength,xlength); %efficient E-L transfer function
    if isnan(mean(U3(:,3)))==1 || isnan(mean(U3(:,2)))==1 ||isnan(mean(U3(:,1)))==1 
       'error, Us not computed'
       break   
    end
    
    U3(stuckb,:)=0;
    U3(stuckt,3)=uz2(Em,2,3); U3(stuckt,1:2)=0;
    %calculate new position from the new velocity using Euler's Method
    Xback=zeros(numOfPoints,1);
    X3=U3*dt/2/st+X;
        
    Xpast(X3(:,3)>zlength)=Xpast(X3(:,3)>zlength)+zlength;
    Xback(X3(:,3)<0)=Xback(X3(:,3)<0)-zlength;
    X3(:,3)=mod(X3(:,3),zlength);
    X3(:,1)=mod(X3(:,1),xlength);
    X3(stuckb,2)=0;
    X3(stuckt,2)=ylength;
    X3((X3(:,2)>ylength),2)=ylength;
    X3((X3(:,2)<0),2)=0;
    
    X3(lower,2)=0;
    X3(upper,2)=ylength;
    
    
    Zdistr=zeros(size(Zdist));
    Zdistl=Zdistr;
    Xdistr=Zdistr; 
    Xdistl=Zdistl;
    %calculate forces at the new positions
    
%     Xcalc=[X(:,1),X(:,2),X(:,3)+zlength*Xpast];
    Xdist(matind)=X2(colind,1)-X3(rowind,1);
    Xdistr(matind)=X2(colind,1)-X3(rowind,1)+xlength;
    Xdistl(matind)=X2(colind,1)-X3(rowind,1)-xlength;
    Ydist(matind)=X2(colind,2)-X3(rowind,2);
    Zdist(matind)=X2(colind,3)-X3(rowind,3);
    Zdistr(matind)=X2(colind,3)-X3(rowind,3)+zlength;
    Zdistl(matind)=X2(colind,3)-X3(rowind,3)-zlength;
    
    zind=find(abs(Zdist)>abs(Zdistr));
    Zdist(zind)=Zdistr(zind);
    zind=find(abs(Zdist)>abs(Zdistl));
    Zdist(zind)=Zdistl(zind);
    
    xind=find(abs(Xdist)>abs(Xdistr));
    Xdist(xind)=Xdistr(xind);
    xind=find(abs(Xdist)>abs(Xdistl));
    Xdist(xind)=Xdistl(xind);
    
    DSq=Xdist.^2+Ydist.^2+Zdist.^2;
    D=sqrt(DSq);
    
    breakind=find(D>2*d0);
    if ~isempty(breakind)

        A(breakind)=0;
        [rowind,colind]=find(A);
        matind=find(A);
        Break=Break+1;
    end
    Dtemp=(D-d0).*K./D;
    Fx=Xdist.*Dtemp;
    Fy=Ydist.*Dtemp;
    Fz=Zdist.*Dtemp;
    
    Fx(~A)=0;
    Fy(~A)=0;
    Fz(~A)=0;
    Fx(isnan(Fx))=0;
    Fy(isnan(Fy))=0;
    Fz(isnan(Fz))=0;
    XrelVel(matind)=v0*(U3(colind,1)-U3(rowind,1));
    YrelVel(matind)=v0*(U3(colind,2)-U3(rowind,2));
    ZrelVel(matind)=v0*(U3(colind,3)-U3(rowind,3));    
    
    VrelDotXrel=b.*(XrelVel.*Xdist+YrelVel.*Ydist+ZrelVel.*Zdist)./DSq;
    XFb=VrelDotXrel.*Xdist;
    YFb=VrelDotXrel.*Ydist;
    ZFb=VrelDotXrel.*Zdist;

%     XFb=b*XrelVel;
%     YFb=b*YrelVel;
%     ZFb=b*ZrelVel;

    XFb(~A)=0;
    YFb(~A)=0;
    ZFb(~A)=0;
    XFb(isnan(XFb))=0;
    YFb(isnan(YFb))=0;
    ZFb(isnan(ZFb))=0;
   
    Fx(stuckt,:)=0;
    Fy(stuckt,:)=0;
    Fz(stuckt,:)=0;
    XFb(stuckb,:)=0;
    YFb(stuckb,:)=0;
    ZFb(stuckb,:)=0;

    Fx(stuckb,:)=0;
    Fy(stuckb,:)=0;
    Fz(stuckb,:)=0;
    XFb(stuckt,:)=0;
    YFb(stuckt,:)=0;
    ZFb(stuckt,:)=0;
    
    %Transfer the new forces and densities to the Eulerian points
    [Efx3,Efy3,Efz3]=transferLtoE3Dper_e2(h,x,y,z,X3,sum(Fx,2)+sum(XFb,2),sum(Fy,2)+sum(YFb,2),sum(Fz,2)+sum(ZFb,2),d0mean,zlength,xlength);


    if addlvisc>0 || addldens>0
        viscmat{1}=transferLtoEvisc3Dper_e2(h,x,y,z,X3,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
        
     Edens{1}.Edensin=transferLtoEdens3Dper_e2(h,x,y,z,X3,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);

        for c5=2:levelsV
            [Em1,En1,Ep1]=size(Edens{c5}.x2h);
            
            
            viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
            viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
            viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
            viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;

            Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
            Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
            Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
            Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;
       


        end
    end
        
        
[vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
    

[uhalfx,uhalfy,uhalfz,verr(c3)]=multigVEL3D_CGper2(dt,h,ux2,uy2,uz2,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx3*fc,Efy3*fc,Efz3*fc,...
    Em,En,Ep,levelsV,vcoef,st,eust,Re);

    if size(p)~=[Em,En,Ep];
        p=zeros(Em,En,Ep);
        issue=1+issue;
    end
     if size(uhalfx)~=[Em,En,Ep];
        uxhalf=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfy)~=[Em,En,Ep];
        uhalfy=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfz)~=[Em,En,Ep];
        uhalfz=zeros(Em,En,Ep);
        issue=1+issue;
     end
        
    %Now solve for new pressure
    [p,perr(c3)]=multigPRESSUREprod3Dper2(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levelsP,pcoef,eust);
   
    %Solve for the new velocity profile u
    
%     parfor 


edenscoef=eust*dt/(2*h)./Edens{1}.Edensin(2:Emm,2:En,2:Ep);
 
    ux3(2:Emm,2:Enm,2:Epm)=(p(2:Emm,1:Enmm,2:Epm)-...
        p(2:Emm,3:En,2:Epm)).*edenscoef(:,1:Enmm,1:Epmm)+uhalfx(2:Emm,2:Enm,2:Epm);
    
    ux3(2:Emm,2:Enm,Ep)=(p(2:Emm,1:Enmm,Ep)-...
        p(2:Emm,3:En,Ep)).*edenscoef(:,1:Enmm,Epm)+uhalfx(2:Emm,2:Enm,Ep);
    ux3(2:Emm,2:Enm,1)=ux3(2:Emm,2:Enm,Ep);
   
    ux3(2:Emm,En,2:Epm)=(p(2:Emm,Enm,2:Epm)-...
        p(2:Emm,2,2:Epm)).*edenscoef(:,Enm,1:Epmm)+uhalfx(2:Emm,En,2:Epm);
    ux3(2:Emm,1,2:Epm)=ux3(2:Emm,En,2:Epm);
    
    ux3(2:Emm,En,Ep)=(p(2:Emm,Enm,Ep)-p(2:Emm,2,Ep)).*edenscoef(:,Enm,Epm)+uhalfx(2:Emm,En,Ep);
    ux3(2:Emm,En,1)=ux3(2:Emm,En,Ep);
    ux3(2:Emm,1,Ep)=ux3(2:Emm,En,Ep);
    ux3(2:Emm,1,1)=ux3(2:Emm,En,Ep);
    
    
    uy3(2:Emm,2:En,2:Epm)=(p(1:Emmm,2:En,2:Epm)-...
        p(3:Em,2:En,2:Epm)).*edenscoef(:,:,1:Epmm)+uhalfy(2:Emm,2:En,2:Epm);
    uy3(2:Emm,2:En,Ep)=(p(1:Emmm,2:En,Ep)-...
        p(3:Em,2:En,Ep)).*edenscoef(:,:,Epm)+uhalfy(2:Emm,2:En,Ep);
    uy3(2:Emm,2:En,1)=uy3(2:Emm,2:En,Ep);
    uy3(2:Emm,1,:)=uy3(2:Emm,En,:);
    
    uz3(2:Emm,2:En,2:Epm)=(p(2:Emm,2:En,1:Epmm)-...
        p(2:Emm,2:En,3:Ep)).*edenscoef(:,:,1:Epmm)+uhalfz(2:Emm,2:En,2:Epm);
    uz3(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Epm)-...
        p(2:Emm,2:En,2)).*edenscoef(:,:,Epm)+uhalfz(2:Emm,2:En,Ep);
    uz3(2:Emm,2:En,1)=uz3(2:Emm,2:En,Ep);
    uz3(2:Emm,1,:)=uz3(2:Emm,En,:);











  U=transferEtoLvel3Dper_e2a(h,ux3,uy3,uz3,x,y,z,X3,zlength,xlength); %efficient E-L transfer function
    if isnan(mean(U2(:,3)))==1 || isnan(mean(U2(:,2)))==1 ||isnan(mean(U2(:,1)))==1 
       'error, U2 not computed'
       break   
    end
    
    U(stuckb,:)=0;
    U(stuckt,3)=uz3(Em,2,3); U(stuckt,1:2)=0;
    %calculate new position from the new velocity using Euler's Method
    Xback=zeros(numOfPoints,1);
        
    X=X+dt/6*(U1+2*U2+2*U3+U);
    Xpast(X(:,3)>zlength)=Xpast(X(:,3)>zlength)+zlength;
    Xback(X(:,3)<0)=Xback(X(:,3)<0)-zlength;
    X(:,3)=mod(X(:,3),zlength);
    X(:,1)=mod(X(:,1),xlength);
    X(stuckb,2)=0;
    X(stuckt,2)=ylength;
    X((X(:,2)>ylength),2)=ylength;
    X((X(:,2)<0),2)=0;
    
    X(lower,2)=0;
    X(upper,2)=ylength;
    
    
    Zdistr=zeros(size(Zdist));
    Zdistl=Zdistr;
    Xdistr=Zdistr; 
    Xdistl=Zdistl;
    %calculate forces at the new positions
    
%     Xcalc=[X(:,1),X(:,2),X(:,3)+zlength*Xpast];
    Xdist(matind)=X(colind,1)-X(rowind,1);
    Xdistr(matind)=X(colind,1)-X(rowind,1)+xlength;
    Xdistl(matind)=X(colind,1)-X(rowind,1)-xlength;
    Ydist(matind)=X(colind,2)-X(rowind,2);
    Zdist(matind)=X(colind,3)-X(rowind,3);
    Zdistr(matind)=X(colind,3)-X(rowind,3)+zlength;
    Zdistl(matind)=X(colind,3)-X(rowind,3)-zlength;
    
    zind=find(abs(Zdist)>abs(Zdistr));
    Zdist(zind)=Zdistr(zind);
    zind=find(abs(Zdist)>abs(Zdistl));
    Zdist(zind)=Zdistl(zind);
    
    xind=find(abs(Xdist)>abs(Xdistr));
    Xdist(xind)=Xdistr(xind);
    xind=find(abs(Xdist)>abs(Xdistl));
    Xdist(xind)=Xdistl(xind);
    
    DSq=Xdist.^2+Ydist.^2+Zdist.^2;
    D=sqrt(DSq);
    
    breakind=find(D>2*d0);
    if ~isempty(breakind)

        A(breakind)=0;
        [rowind,colind]=find(A);
        matind=find(A);
        Break=Break+1;
    end
    Dtemp=(D-d0).*K./D;
    Fx=Xdist.*Dtemp;
    Fy=Ydist.*Dtemp;
    Fz=Zdist.*Dtemp;
    
    Fx(~A)=0;
    Fy(~A)=0;
    Fz(~A)=0;
    Fx(isnan(Fx))=0;
    Fy(isnan(Fy))=0;
    Fz(isnan(Fz))=0;
    XrelVel(matind)=v0*(U(colind,1)-U(rowind,1));
    YrelVel(matind)=v0*(U(colind,2)-U(rowind,2));
    ZrelVel(matind)=v0*(U(colind,3)-U(rowind,3));    
    
    VrelDotXrel=b.*(XrelVel.*Xdist+YrelVel.*Ydist+ZrelVel.*Zdist)./DSq;
    XFb=VrelDotXrel.*Xdist;
    YFb=VrelDotXrel.*Ydist;
    ZFb=VrelDotXrel.*Zdist;

%     XFb=b*XrelVel;
%     YFb=b*YrelVel;
%     ZFb=b*ZrelVel;

    XFb(~A)=0;
    YFb(~A)=0;
    ZFb(~A)=0;
    XFb(isnan(XFb))=0;
    YFb(isnan(YFb))=0;
    ZFb(isnan(ZFb))=0;
	
	
    [Efxt,Efyt,Efzt]=transferLtoE3Dper_e2(h,x,y,z,X,sum(Fx,2)+sum(XFb,2),sum(Fy,2)+sum(YFb,2),sum(Fz,2)+sum(ZFb,2),d0mean,zlength,xlength);


    Fx(stuckt,:)=0;
    Fy(stuckt,:)=0;
    Fz(stuckt,:)=0;
    XFb(stuckb,:)=0;
    YFb(stuckb,:)=0;
    ZFb(stuckb,:)=0;

    Fx(stuckb,:)=0;
    Fy(stuckb,:)=0;
    Fz(stuckb,:)=0;
    XFb(stuckt,:)=0;
    YFb(stuckt,:)=0;
    ZFb(stuckt,:)=0;
    
    %Transfer the new forces and densities to the Eulerian points
    [Efx,Efy,Efz]=transferLtoE3Dper_e2(h,x,y,z,X,sum(Fx,2)+sum(XFb,2),sum(Fy,2)+sum(YFb,2),sum(Fz,2)+sum(ZFb,2),d0mean,zlength,xlength);


    if addlvisc>0 || addldens>0
        viscmat{1}=transferLtoEvisc3Dper_e2(h,x,y,z,X,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
        
     Edens{1}.Edensin=transferLtoEdens3Dper_e2(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);

        for c5=2:levelsV
            [Em1,En1,Ep1]=size(Edens{c5}.x2h);
            
            
            viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
            viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
            viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
            viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;

            Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
            Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
            Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
            Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;
       


        end
    end
        
        
[vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
    

[uhalfx,uhalfy,uhalfz,verr(c3)]=multigVEL3D_CGper2(dt,h,ux,uy,uz,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx*fc,Efy*fc,Efz*fc,...
    Em,En,Ep,levelsV,vcoef,st,eust,Re);

    if size(p)~=[Em,En,Ep];
        p=zeros(Em,En,Ep);
        issue=1+issue;
    end
     if size(uhalfx)~=[Em,En,Ep];
        uxhalf=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfy)~=[Em,En,Ep];
        uhalfy=zeros(Em,En,Ep);
        issue=1+issue;
     end
     if size(uhalfz)~=[Em,En,Ep];
        uhalfz=zeros(Em,En,Ep);
        issue=1+issue;
     end
        
    %Now solve for new pressure
    [p,perr(c3)]=multigPRESSUREprod3Dper2(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levelsP,pcoef,eust);
   
    %Solve for the new velocity profile u
    
%     parfor 


edenscoef=eust*dt/(2*h)./Edens{1}.Edensin(2:Emm,2:En,2:Ep);
 
    ux(2:Emm,2:Enm,2:Epm)=(p(2:Emm,1:Enmm,2:Epm)-...
        p(2:Emm,3:En,2:Epm)).*edenscoef(:,1:Enmm,1:Epmm)+uhalfx(2:Emm,2:Enm,2:Epm);
    
    ux(2:Emm,2:Enm,Ep)=(p(2:Emm,1:Enmm,Ep)-...
        p(2:Emm,3:En,Ep)).*edenscoef(:,1:Enmm,Epm)+uhalfx(2:Emm,2:Enm,Ep);
    ux(2:Emm,2:Enm,1)=ux(2:Emm,2:Enm,Ep);
   
    ux(2:Emm,En,2:Epm)=(p(2:Emm,Enm,2:Epm)-...
        p(2:Emm,2,2:Epm)).*edenscoef(:,Enm,1:Epmm)+uhalfx(2:Emm,En,2:Epm);
    ux(2:Emm,1,2:Epm)=ux(2:Emm,En,2:Epm);
    
    ux(2:Emm,En,Ep)=(p(2:Emm,Enm,Ep)-p(2:Emm,2,Ep)).*edenscoef(:,Enm,Epm)+uhalfx(2:Emm,En,Ep);
    ux(2:Emm,En,1)=ux(2:Emm,En,Ep);
    ux(2:Emm,1,Ep)=ux(2:Emm,En,Ep);
    ux(2:Emm,1,1)=ux(2:Emm,En,Ep);
    
    
    uy(2:Emm,2:En,2:Epm)=(p(1:Emmm,2:En,2:Epm)-...
        p(3:Em,2:En,2:Epm)).*edenscoef(:,:,1:Epmm)+uhalfy(2:Emm,2:En,2:Epm);
    uy(2:Emm,2:En,Ep)=(p(1:Emmm,2:En,Ep)-...
        p(3:Em,2:En,Ep)).*edenscoef(:,:,Epm)+uhalfy(2:Emm,2:En,Ep);
    uy(2:Emm,2:En,1)=uy(2:Emm,2:En,Ep);
    uy(2:Emm,1,:)=uy(2:Emm,En,:);
    
    uz(2:Emm,2:En,2:Epm)=(p(2:Emm,2:En,1:Epmm)-...
        p(2:Emm,2:En,3:Ep)).*edenscoef(:,:,1:Epmm)+uhalfz(2:Emm,2:En,2:Epm);
    uz(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Epm)-...
        p(2:Emm,2:En,2)).*edenscoef(:,:,Epm)+uhalfz(2:Emm,2:En,Ep);
    uz(2:Emm,2:En,1)=uz(2:Emm,2:En,Ep);
    uz(2:Emm,1,:)=uz(2:Emm,En,:);



    %Calculate Stresses on top wall.   
    eShear(c3) = sum(sum(sum(Efzt(:,1:Enm,1:Epm))))*dx^3/(0.81)*charLength;  %S.n = t=F/A --- *charLength??
    vShear(c3) =( sum(sum(viscmat{1}(Em,1:Enm,1:Epm).*(3*uz(Em,1:Enm,1:Epm)-4*uz(Emm,1:Enm,1:Epm)+uz(Emmm,1:Enm,1:Epm))))/(2*dx)...
        +sum(sum(sum(viscmat{1}(Em-10:Emm,1:Enm,1:Epm).*(uz(Em-9:Em,1:Enm,1:Epm)-uz(Em-11:Emmm,1:Enm,1:Epm)))))/(2*dx))*dx^3/(xlength*10*h*zlength);
    Shear_Force(c3)=-eShear(c3)*dx^3-vShear(c3)*dx^3;

    [bStrain(c3), fStrain(c3),S,EU,bStraintop(c3),fStraintop(c3),Shear_Strainb,Shear_Strainf] = Calculate_Strains2(A,X,ux,uy,uz,dx,X0,S,dt,st,xlength,ylength,zlength,upper,h,x,y,z,numOfnonzero);

    tStrain(c3)=bStrain(c3)+fStrain(c3);
    eShear(c3);
    v0*visc0/charLength*vShear(c3);

    
%     profile off
%     profile viewer
    toc(tstart)
    if mod(c3,10000)==0
        uxsave=[uxsave, ux];
        uysave=[uysave,uy];
        uzsave=[uzsave,uz];
        psave=[psave, p];
        Xstore=[Xstore, X];
        
%         vsave=[vsave, viscmat{1}];
%         densave=[densave, Edens{1}.Edensin]; 
%         Ssave{1}=[Ssave{1},S{1}];
%         Ssave{2}=[Ssave{2},S{2}];
%         Ssave{3}=[Ssave{3},S{3}];
%         save sim3DShroom0_0dens500_0viscFromDataSTW3fmax4dt0001dx1_256_dt_0001_t_01_scl_rc_1.mat
    end
%    if mod(c3,20)==0
%        plot(vShear(1:c3))
%        pause(1)
%    end

    if mod(c3,200)==0
    	vShear(c3)
	eShear(c3)
        str1=num2str(w);
	str2=num2str(fmax/1000);
	str3=num2str(B);
	str4=num2str(1/dt);
	str5=num2str(c3/100);
        runid=['Sim_',str1,'_',str2,'_',str3,'_',str4];
        save(['/home/jstotsky/scratch/',runid,'.mat']);
    end
end

% matlabpool close
% save ShroomFromData.mat
