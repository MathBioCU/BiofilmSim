levelsV=5;
levelsP=5;


xlength=0.9397*2;%diameter of the tube
ylength=0.9397*2;%radius of the tube
zlength=0.9397;%length of the tube for the computational domain

%%%%%%%%%%%%%dimensional constants
% w=1;

co = min(4, abs(0.09/(0.0108*log(w)+0.0465)));%correction to keep strain amp const over various frequencies
e0=ylength/1.8*4*4.5*10*10^-6*tan(0.13)*co;
v0=e0*w;%speed at the middle
visc0=10^-3;%this is dynamic viscosity of water, units of kg/m/s
rho0=998;%9.983*10^(-7);%in units of kg/m^3
charLength=10*10^-6;
timefreq=1;
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

h=0.9397/64;%64 should be min, this is the hx hy, and hz
% n=1/h-1;
% levels=5;
% dt=.01;
xvec=0:h:xlength;
yvec=0:h:ylength;
zvec=0:h:zlength;

[x,y,z]=meshgrid(xvec,yvec,zvec);

tend=numtimesteps*dt;
% connectdist=3/18;

[Xb,Yb,Zb]=importfile('37 Deg Test 1, Live Cells.txt');
%x1=[Xb, Yb,Zb];
%x1=x1(x1(:,1)<30 & x1(:,2)<30 & x1(:,3)<15 & x1(:,2)>0 & x1(:,1)>0 ,:);
%x1=[x1(:,1)-5,x1(:,2)-5,x1(:,3)]; %set edge from 10 microns to 0 microns
%X=[x1(:,1),x1(:,2),x1(:,3)]*mult/charLength;
X=[Xb,Zb,Yb]*mult/charLength;
bthickness=0.07; %thickness of region where bacteria are adhered to walls
% X=X(find(X(:,2)<0.2),:);
X=X(X(:,1)>0,:);
X=X(X(:,1)<xlength,:);
X=X(X(:,2)>0,:);
X=X(X(:,2)<zlength,:);
X(X(:,3)>ylength & X(:,3)<ylength+bthickness/2,3)=ylength;
X(X(:,3)<0 & X(:,3)>-bthickness/3,3)=0;
X=X(X(:,3)<=ylength,:);
X=X(X(:,3)>=0,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to perform rolling tests, change allowed range for X(:,3)
% to something like X(:,3)>0.7 and X(:,3)<1.3
%so that all of the biofilm nodes in are center of flow
%may also want to restrict x and z directions too

%X=X(X(:,1)<xlength-0.2,:);
%X=X(X(:,1)>0.2,:);
%X=X(X(:,2)>0.2,:);
%X=X(X(:,2)<zlength-0.2,:);
%X=X(X(:,3)<=ylength-0.7,:);
%X=X(X(:,3)>=0.7,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X=X(X(:,2)>0.01,:);
Xupper=X(X(:,3)>=ylength-bthickness,:);
%Xupper(:,3)=ylength;
%X=X(X(:,3)<ylength,:);
%X=[X;Xupper];

Xlower=X(X(:,3)<=bthickness+0.03,:);
%Xlower(:,3)=0;
%X=X(X(:,3)>=0,:);
%X=[X;Xlower];


% X=[0,0,0];


lower=find(X(:,3)<=bthickness+0.03);
upper=find(X(:,3)>=ylength-bthickness);

% load X.mat
% load X1.mat
[X1,Y1,Z1]=meshgrid(0:d0mean:0.2,0:d0mean:0.2,0:d0mean:0.2);
Xt=[X1(:),Y1(:),Z1(:)];
Xt=Xt*charLength;

X=[X(:,1),X(:,3),X(:,2)];



numOfPoints=size(X,1);
numOfnonzero=numOfPoints;%-length(find(X(:,3)<=0.05))-length(find(X(:,3)>=ylength-0.05));
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
Xrv=Xdist; Yrv=Xdist; Zrv=Xdist;
[Em,En,Ep]=size(x);


XFe=zeros(size(A));
YFe=zeros(size(A));
ZFe=zeros(size(A));
    
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
E=E./d0;

E(~A)=0;
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

crossbarrier=zeros(numtimesteps,1);
[vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
stuckb=find(X(:,2)<bthickness+0.03);
stuckt=find(X(:,2)>ylength-bthickness);
Xpast=zeros(numOfPoints,1);
%% Go through all of the timesteps
%lower=[];
%stuckb=[];

uxsave=ux;
uysave=uy;
uzsave=uz;
psave=p;
vsave=viscmat{1};
densave=Edens{1}.Edensin;

issue=1;

S=zeros(108,3);
for i=1:6
    for j=1:6
        S(6*(i-1)+j,:)=[i*xlength/7, ylength-2*h, j*zlength/7];
    end
end

S(37:72,:)=S(1:36,:);
S(73:108,:)=S(1:36,:);
S(37:72,2)=S(37:72,2)-h;
S(73:108,2)=S(73:108,2)-2*h;

Sstore=S;
S0=S;
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

dx=h;
    
uz(1,:,:)=0;
uz(Em,:,:)=(exp(2*c3*w*dt)-1)*((exp(4*c3*w*dt)-1)*cos(c3*dt)+exp(2*c3*w*dt)*8*sin(c3*w*dt))/(4*(1+exp(2*c3*w*dt))^3);
uzs(1,:,:)=uz(1,:,:);
uzs(Em,:,:)=uz(Em,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('~/scratch/w62.83_f380_b0_cnd185_u500_dt500Profile2.mat');
%numtimesteps=10+10*numtimesteps;
%dt=dt/w/10;
%timefreq=1;
%st=timefreq*charLength/v0;
%eust=eu/st;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c3=1:numtimesteps   %change to 1 for actual simulation
    
    c3
  [vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);


    uxm=ux;
    uym=uy;
    uzm=uz;
    
    dx=h;
    tstart=tic;
   
    
    uz(1,:,:)=0;
    uz(Em,:,:)=(exp(2*c3*w*dt)-1)*((exp(4*c3*w*dt)-1)*cos(c3*w*dt)+exp(2*c3*w*dt)*8*sin(c3*w*dt))/(4*(1+exp(2*c3*w*dt))^3);
    uzs(1,:,:)=uz(1,:,:);
    uzs(Em,:,:)=uz(Em,:,:);

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
    

    %Transfer the velocity to the Lagrangian points
    Um=U;
    U=transferEtoLvel3Dper_e2a(h,ux,uy,uz,x,y,z,X,zlength,xlength); %efficient E-L transfer function
    if isnan(mean(U(:,3)))==1 || isnan(mean(U(:,2)))==1 ||isnan(mean(U(:,1)))==1 
       'error, U not computed'
       break   
    end
    
    U(stuckb,3)=0;
    U(stuckt,3)=uz(Em,2,3); %U(stuckt,1:2)=0;
    %calculate new position from the new velocity using Euler's Method
    Xback=zeros(numOfPoints,1);
    X=U*dt/st+X;

    %lower=find(X(:,2)<=bthickness);
    %upper=find(X(:,2)>=ylength-bthickness);
    %stuckb=find(X(:,2)<=bthickness);
    %stuckt=find(X(:,2)>=ylength-bthickness);

    Xpast(X(:,3)>zlength)=Xpast(X(:,3)>zlength)+zlength;
    Xback(X(:,3)<0)=Xback(X(:,3)<0)-zlength;
    X(:,3)=mod(X(:,3),zlength);
    X(:,1)=mod(X(:,1),xlength);
%    X(stuckb,2)=0;
%    X(stuckt,2)=ylength;
    X((X(:,2)>=ylength),2)=ylength;
    X((X(:,2)<=0),2)=0;
    
%    X(lower,2)=0;
%    X(upper,2)=ylength;
    
    
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

    
%    XrelVel(matind)=(abs(v0*(U(colind,1)-U(rowind,1)))).^(0.9).*sign(U(colind,1)-U(rowind,1));
%    YrelVel(matind)=(abs(v0*(U(colind,2)-U(rowind,2)))).^(0.9).*sign(U(colind,2)-U(rowind,2));
%    ZrelVel(matind)=(abs(v0*(U(colind,3)-U(rowind,3)))).^(0.9).*sign(U(colind,3)-U(rowind,3));    
  
    Xrv(matind)=v0*(U(colind,1)-U(rowind,1));
    Yrv(matind)=v0*(U(colind,2)-U(rowind,2));
    Zrv(matind)=v0*(U(colind,3)-U(rowind,3));    
    aa=0.1;
    XrelVel=(1-(1./(1+abs(Xrv).^(1-aa))).*exp(-abs(Xrv).^(1-aa))).*abs(Xrv).^aa.*sign(Xrv);
    YrelVel=(1-(1./(1+abs(Yrv).^(1-aa))).*exp(-abs(Yrv).^(1-aa))).*abs(Yrv).^aa.*sign(Yrv);
    ZrelVel=(1-(1./(1+abs(Zrv).^(1-aa))).*exp(-abs(Zrv).^(1-aa))).*abs(Zrv).^aa.*sign(Zrv);

    VrelDotXrel=(XrelVel.*Xdist+YrelVel.*Ydist+ZrelVel.*Zdist)./DSq;
    XFb=b.*(VrelDotXrel.*Xdist);
    YFb=b.*(VrelDotXrel.*Ydist);
    ZFb=b.*(VrelDotXrel.*Zdist);


%     XFb=b*XrelVel;
%     YFb=b*YrelVel;
%     ZFb=b*ZrelVel;

    XFb(~A)=0;
    YFb(~A)=0;
    ZFb(~A)=0;
    XFb(isnan(XFb))=0;
    YFb(isnan(YFb))=0;
    ZFb(isnan(ZFb))=0;

    XFe=1./(1+E/dt).*(Fx+XFb+E/dt.*XFe);
    YFe=1./(1+E/dt).*(Fy+YFb+E/dt.*YFe);
    ZFe=1./(1+E/dt).*(Fz+ZFb+E/dt.*ZFe);

    XFe(~A)=0;
    YFe(~A)=0;
    ZFe(~A)=0;
    XFe(isnan(XFe))=0;
    YFe(isnan(YFe))=0;
    ZFe(isnan(ZFe))=0;
  
    [Efxt,Efyt,Efzt]=transferLtoE3Dper_e2a(h,x,y,z,X(stuckt,:),sum(XFe(stuckt,:),2),sum(YFe(stuckt,:),2),sum(ZFe(stuckt,:),2),d0mean,zlength,xlength);

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
    [Efx,Efy,Efz]=transferLtoE3Dper_e2a(h,x,y,z,X,sum(XFe,2),sum(YFe,2),sum(ZFe,2),d0mean,zlength,xlength);


    if addlvisc>0 || addldens>0
        viscmat{1}=transferLtoEvisc3Dper_e2a(h,x,y,z,X,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
        
     Edens{1}.Edensin=transferLtoEdens3Dper_e2a(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);

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
    

    %Calculate Stresses on top wall.   
    eShear(c3) = sum(sum(sum(Efzt(:,1:Enm,1:Epm))))*dx^3/(0.81)*charLength;  %S.n = t=F/A --- *charLength??
    vShear(c3) =sum(sum(sum(viscmat{1}(Em-10:Emm,1:Enm,1:Epm).*(uz(Em-9:Em,1:Enm,1:Epm)-uz(Em-11:Emmm,1:Enm,1:Epm)))))/(2*dx)*dx^3/(xlength*10*h*zlength);
    Shear_Force(c3)=-eShear(c3)*dx^3-vShear(c3)*dx^3;

    [bStrain(c3), fStrain(c3),S] = Calculate_Strains2(A,X,S,ux,uy,uz,dx,X0,S0,dt,st,xlength,ylength,zlength,upper,h,x,y,z,numOfnonzero);

    tStrain(c3)=bStrain(c3)+fStrain(c3);
    eShear(c3);
    v0*visc0/charLength*vShear(c3);

%     profile off
%     profile viewer
    Sstore=[Sstore,S];
    Xstore=[Xstore,X];
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
    if mod(c3,100)==0
        str1=num2str(w);
	str2=num2str(fmax/1000);
	str3=num2str(B01);
	str4=num2str(1/(dt*w));
	str5=num2str(c3/100);
        runid=['Sim_',str1,'_',str2,'_',str3,'_',str4];
        save(['/lustre/janus_scratch/jast0817/Simulation_Test_Parameters/',runid,'.mat']);
    end
    	
end

% matlabpool close
% save ShroomFromData.mat
