%Condense portions of code used to updating velocity and position, work on
%different time stepping methods


levelsV=5;
levelsP=5;

% B01=0;
% b=0;
% E=0;
% w=62.83;
% dt=1/(500*w);
% numtimesteps=100;
% fmax=25000;
% connectdist=0.1622;
% addlvisc=250;

mult=10^-6;
charLength=10*10^-6;
[Xb,Yb,Zb]=importfile('~/Documents/Biofilm_Research/Data/37DegTest2LiveCells.txt');
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
e0=ylength/1.8*4*4.5*10*10^-6*tan(0.13)*co/0.90;  
v0=e0*w;%speed at the middle
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

xvec=0:h:xlength;
yvec=0:h:ylength;
zvec=0:h:zlength;

[x,y,z]=meshgrid(xvec,yvec,zvec);

tend=numtimesteps*dt;

bthickness=0.04; %thickness of region where bacteria are adhered to walls

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
    
ux=zeros(Em,En,Ep);
uy=ux;
uz=ux;
uxs=ux; uys=uy; uzs=uz;

%for comparison without biofilm to exact oscillatory solution
 %kw=(w/(2*visc0/rho0))^(1/2);
 %for i=1:Em
 %    uz(i,:,:)=sqrt(real(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2+imag(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2).*sin(angle(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))));
 %end


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
Xrstore=X;
Xreal=X;
X0=X;
%Transfer the density to the Eulerian grid from the Lagrangian/compute
%new densities

[Edens,viscmat, viscmatmid]=dens_visc(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,ylength, xlength,visc,addlvisc,levelsV,Em,En,Ep);
coefmult=1;

crossbarrier=zeros(numtimesteps,1);
[vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
stuckb=find(X(:,2)<bthickness);
stuckt=find(X(:,2)>ylength-bthickness);
Xpast=zeros(numOfPoints,1);
%% Go through all of the timesteps

uxsave=ux;
uysave=uy;
uzsave=uz;
psave=p;
vsave=viscmat{1};
densave=Edens{1}.Edensin;
issue=1;

S=zeros(336,3);
for i=1:10
    for j=1:10
        S(10*(i-1)+j,:)=[i*xlength/11, ylength-2*h, j*zlength/11];
    end
end

S(101:200,:)=S(1:100,:);
S(201:300,:)=S(1:100,:);
S(101:200,2)=S(101:200,2)-h;
S(201:300,2)=S(201:300,2)-2*h;
S(301:336,2)=[0:ylength/35:ylength]';
S(301:336,1)=xlength/2;
S(301:336,3)=zlength/2;


Sstore=S;
S0=S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eShear=zeros(numtimesteps,1); vShear=eShear; Shear_Force=eShear;
perr=zeros(numtimesteps,1); verr=perr; vShear2=eShear;


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
uxm=ux;
uym=uy;
uzm=uz;

dx=h;
    
uz(1,:,:)=0;
uz(Em,:,:)=(exp(2*c3*w*dt)-1)*((exp(4*c3*w*dt)-1)*cos(c3*dt)+exp(2*c3*w*dt)*8*sin(c3*w*dt))/(4*(1+exp(2*c3*w*dt))^3);
uzs(1,:,:)=uz(1,:,:);
uzs(Em,:,:)=uz(Em,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('Sim_1.254_95_0_500.mat'); i.e. to loop data from previous simulation
% c31=c3;
c31=1;

%numtimesteps=10+10*numtimesteps;
%dt=dt/w/10;
%timefreq=1;
%st=timefreq*charLength/v0;
%eust=eu/st;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
time=0;
times=zeros(numtimesteps+1,1);
for i=2:numtimesteps+1
	times(i)=times(i-1)+dts(i-1);
end

speed=zeros(numtimesteps,1);


for c3=1:numtimesteps  
    c3
   %dt=dts(c3);
   time=time+dt;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uxm=ux; uym=uy; uzm=uz;
    dx=h;
    tstart=tic;
    uz(1,:,:)=0;
	%Periodic STrain
    uz(Em,:,:)=(exp(2*time*w)-1)*((exp(4*time*w)-1)*cos(time*w)+exp(2*time*w)*8*sin(time*w))/(4*(1+exp(2*time*w))^3); %sin(w*c3*dt);
    uzs(1,:,:)=uz(1,:,:);
    uzs(Em,:,:)=uz(Em,:,:);

%    Compliance
%    if c3>1
%    	uz(Em,:,:)=uzm(Em,:,:)+dts(c3)*((2./(1+exp(-55*w*time))-1)*sigma0+eShear(c3-1)-v0*visc0/charLength*vShear(c3-1))/(rho0*sum(sum(Edens{1}.Edensin(Em,:,:))))/(charLength*h)/v0; %jumps in dt cause wiggling due to force overshooting
%	%try to adjust dt at jumps in dt to avoid wiggling around there. Done by assuming that dt=dts(c3-1) so that no jump in acceleration suddenly occurs
%    else
%	uz(Em,:,:)=uzm(Em,:,:)+dt*(2./(1+exp(-55*w*time))-1)*sigma0/rho0/charLength/1000;
%   end
    speed(c3)=uz(Em,1,1);
    [vcoef,pcoef]=compute_operator(Edens,levelsV,Re,st,viscmat,viscmatmid,dt);
    [uhalfx,uhalfy,uhalfz,verr(c3)]=multigVEL3D_CGper2(dt,h,ux,uy,uz,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx*fc,Efy*fc,Efz*fc,...
    Em,En,Ep,levelsV,vcoef,st,eust,Re);
        
    %Now solve for new pressure
    [p,perr(c3)]=multigPRESSUREprod3Dper4(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levelsP,pcoef,eust);
    %Solve for the new velocity profile u
    [ux,uy,uz]=DS_Vel(ux,uy,uz,p,eust,dt,h,Edens,Em,En,Ep,uhalfx,uhalfy,uhalfz);
    

    %Transfer the velocity to the Lagrangian points
    Um=U;
    U=transferEtoLvel3Dper_e2a(h,ux,uy,uz,x,y,z,X,zlength,xlength); %efficient E-L transfer function
    if isnan(mean(U(:,3)))==1 || isnan(mean(U(:,2)))==1 ||isnan(mean(U(:,1)))==1 
       fprintf('Error, U not computed\n');
       break   
    end
    
    U(stuckb,3)=0;
    U(stuckt,3)=uz(Em,2,3); %U(stuckt,1:2)=0;
    %calculate new position from the new velocity using Euler's Method
    Xback=zeros(numOfPoints,1);
    X=U*dt/st+X;
    Xreal=U*dt/st+Xreal;
    X(:,3)=mod(X(:,3),zlength);
    X(:,1)=mod(X(:,1),xlength);
    X((X(:,2)>=ylength),2)=ylength;
    X((X(:,2)<=0),2)=0;   

    [XFe,YFe,ZFe,XFet,YFet,ZFet,Break,rowind,colind,matind,Intx,Inty,Intz,Dtemp]=...
       LagrForce(X,xlength,ylength,zlength,matind,colind,rowind,Break,v0,U,Xrv,Yrv,Zrv,dt,K,d0,b,Xdist,Ydist,Zdist,A,XFe,YFe,ZFe,stuckt,stuckb,t(c3),Intx,Inty,Intz,Dtemp);
 
    %Transfer the new forces and densities to the Eulerian points
    [Efx,Efy,Efz]=transferLtoE3Dper_e2a(h,x,y,z,X,sum(XFe,2),sum(YFe,2),sum(ZFe,2),d0mean,zlength,xlength);
    [Efxt,Efyt,Efzt]=transferLtoE3Dper_e2a(h,x,y,z,X(stuckt,:),sum(XFet(stuckt,:),2),sum(YFet(stuckt,:),2),sum(ZFet(stuckt,:),2),d0mean,zlength,xlength);

    [Edens,viscmat,viscmatmid]=DensVisc(addlvisc,addldens,h,x,y,z,X,visc,numOfnonzero,xlength,zlength,Em,En,Ep,d0mean,initdensity,levelsV,viscmat,Edens);        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate Stresses on top wall.   
    eShear(c3) = sum(sum(sum(Efzt(:,1:Enm,1:Epm))))*dx^3/(xlength*zlength)*charLength;  %S.n = t=F/A --- *charLength??
    vShear(c3) =sum(sum(sum(viscmat{1}(Em-11:Emm,1:Enm,1:Epm).*(uz(Em-10:Em,1:Enm,1:Epm)-uz(Em-12:Emmm,1:Enm,1:Epm)))))/(2*dx)*dx^3/(xlength*10*h*zlength);
    vShear2(c3) =sum(sum(sum(viscmat{1}(Em-5:Emm,1:Enm,1:Epm).*(uz(Em-4:Em,1:Enm,1:Epm)-uz(Em-6:Emmm,1:Enm,1:Epm)))))/(2*dx)*dx^3/(xlength*5*h*zlength);

	Shear_Force(c3)=-eShear(c3)*dx^3-vShear(c3)*dx^3;

    [fStrain(c3),S,S0] = Calculate_Strains3(S,ux,uy,uz,S0,dt,st,xlength,ylength,zlength,h,x,y,z);

%    tStrain(c3)=bStrain(c3)+fStrain(c3);
    eShear(c3);
    v0*visc0/charLength*vShear(c3);

%     profile off
%     profile viewer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sstore=[Sstore,S];
    Xstore=[Xstore,X];
    Xrstore=[Xrstore,Xreal];
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
        runid=['testSim_',str1,'_',str2,'_',str3,'_',str4];
        save([runid,'.mat']);
    end
    	
end