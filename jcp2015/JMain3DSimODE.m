%ODE solver for updating biofilm positions, ignores presence of fluid (i.e.
% for very slow fluids where fluid stresses are negligible compared to
% viscoelastic bonds between bacteria) 

mult=10^-6;
charLength=10*10^-6;
[Xb,Yb,Zb]=importfile('37 Deg Test 2, Live Cells.txt');
X=[Xb,Zb,Yb]*mult/charLength;

mdim=min([max(X(:,1))-min(X(:,1)),max(X(:,3))-min(X(:,3)), max(X(:,2))-min(X(:,2))]);

xlength=mdim*1;%width of the tube
ylength=mdim*3;%height of the tube
zlength=mdim;%length of the tube for the computational domain

X(:,3)=X(:,3)-min(X(:,3));
X(:,1)=X(:,1)-min(X(:,1));
X(:,2)=X(:,2)-min(X(:,2));


levelsV=-log2(dx);
levelsP=-log2(dx);
h=mdim*dx;%64 should be min, this is the hx hy, and hz
%%%%%%%%%%%%%dimensional constants

co = min(1.2, abs(0.09/(0.0120*log(w)+0.0465)));% approximate correction to try to keep strain amp const over various frequencies
e0=ylength/1.8*4*4.5*10*10^-6*tan(0.13)*co/1.5; %0.9  
v0=e0*w;%speed at the middle
if ShearRotation==1
    v0=0.001;
end
visc0=10^-3;%this is dynamic viscosity of water, units of kg/m/s
rho0=998;%9.983*10^(-7);%in units of kg/m^3
%v0=sigma0*dt/charLength/h/rho0;
timefreq=1;
d0mean=1.5921*10^-6;%calculated from given Stewart data as (Vol/(# Lagr nodes interior to the Vol))^(1/3)
d0mean=d0mean/charLength;%nondimensionalized to use in the simulations

%%%%%%%%%%%%%%%%%%%%%%
p0=1;
%%%%%%%%%%%%nondimensional constants
st=timefreq*charLength/v0;
eu=p0/rho0/v0^2;
Re=rho0*charLength*v0/visc0;
fc=charLength/rho0/v0^2;
B=b;
eust=eu/st;


%%%%%%%%%%%%%%%%%%%%%%%%

initdensity=1;
addldens=0.12;%
visc=1;

tend=numtimesteps*dt;

bthickness=0.04; %thickness of region where bacteria are adhered to walls

X=X(X(:,1)>0.2,:);
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

if ShearRotation==1
    zlength=zlength*2;
    xvec=0:h:xlength;
    yvec=0:h:ylength;
    zvec=0:h:zlength;

    X(:,2)=X(:,2)+0.4;
    X=X(X(:,1)<xlength-0.2,:);
    X=X(X(:,1)>0.2,:);
    X=X(X(:,2)>0.4,:);
    X=X(X(:,2)<zlength-0.8,:);
    X=X(X(:,3)<=ylength-0.4,:);
    X=X(X(:,3)>=0.4,:);
    mX=mean(X);
    Z=X-ones(length(X),1)*mX;
    Z=Z*pca(Z);
    Z=Z+ones(length(X),1)*[xlength/2,zlength/2,ylength/2];
    X=Z;
    X=X(1:110,:);
    
   
    X=X(X(:,1)<xlength-0.2,:);
    X=X(X(:,1)>0.2,:);
    X=X(X(:,2)>0.4,:);
    X=X(X(:,2)<zlength-0.8,:);
    X=X(X(:,3)<=ylength-0.4,:);
    X=X(X(:,3)>=0.4,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Xupper=X(X(:,3)>=ylength-bthickness,:);

Xlower=X(X(:,3)<=bthickness,:);

%X=[0,0,0];% use to test without biofilm


lower=find(X(:,3)<=bthickness);
upper=find(X(:,3)>=ylength-bthickness);

X=[X(:,1),X(:,3),X(:,2)]; %highlight to get back to original z-orientation
%X=[X(:,2),X(:,3),X(:,1)]; %changed so that velocity is in X orientation


numOfPoints=size(X,1);
numOfnonzero=numOfPoints;
%% initialize all the variables
%A is connectivity matrix for Lagrangian points
A=zeros(numOfPoints);
d0=A;
% D=zeros(numOfPoints);
Xdist=A;
Ydist=A;
Zdist=A;
Intx=A*0;
Inty=A*0;
Intz=A*0;
Dtemp=A*0;
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
XFet=XFe;
YFet=YFe;
ZFet=ZFe;

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
Xrstore=X;
Xreal=X;
X0=X;
%Transfer the density to the Eulerian grid from the Lagrangian/compute
%new densities


crossbarrier=zeros(numtimesteps,1);
stuckb=find(X(:,2)<bthickness);
stuckt=find(X(:,2)>ylength-bthickness);
Xpast=zeros(numOfPoints,1);
%% Go through all of the timesteps

issue=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eShear=zeros(numtimesteps,1); vShear=eShear; Shear_Force=eShear;

t_star=floor(pi/w/dt);
e0=e0/charLength;

Strain=e0*sin([0:numtimesteps]*dt*w).*(1./(1+exp(-2*[0:numtimesteps]*dt*w))-1/2).^2;

fStrain=zeros(numtimesteps,1);
bStrain=fStrain;
tStrain=fStrain;
bStraintop=bStrain;
fStraintop=fStrain;

Break=0;
c3=1;

c31=1;
%for compliance test, helpful to change dt to capture short and long
%dynamics
dts=zeros(numtimesteps,1);
for i=1:180
	dts(i)=dt/10;
end
for i=181:200
	dts(i)=dt/(10-(i-180)/4);
end
dts(201:300)=dt/5;
for i=301:340
	dts(i)=dt/(5-(i-300)/10);
end
dts(341:end)=dt;

for i=1:numtimesteps
    dts(i)=dt/100*i;
end

time=0;
times=zeros(numtimesteps,1);
for i=2:numtimesteps
	times(i)=times(i-1)+dts(i-1);
end
speed=zeros(numtimesteps,1);


for c3=1:numtimesteps
    [XFe,YFe,ZFe,XFet,YFet,ZFet,Break,rowind,colind,matind,Intx,Inty,Intz,Dtemp]=...
        LagrForce(X,xlength,ylength,zlength,matind,colind,rowind,Break,v0,U,...
        Xrv,Yrv,Zrv,dt,K,d0,b,Xdist,Ydist,Zdist,A,XFe,YFe,ZFe,stuckt,stuckb,t,Intx,Inty,Intz,Dtempm);
    
    X=X+[XFe,YFe,ZFe]/m*dt;
end
    
    
    
    
