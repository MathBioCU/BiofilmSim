function [vxsave,vysave,vzsave,vzerr]=multigVEL3D_CGper2(dt,h,ux,uy,uz,uxm,uym,uzm,uhalfx,uhalfy,uhalfz,p,Edens,Efx,Efy,Efz,...
    Em,En,Ep,levels,vcoef,st,eust,Re)
allowerr=5*10^-9;
% allowerr=10^-11;
numofc=inf;
numofFMV=5;
% En=En+1;
Emm=Em-1;
Enm=En-1;
Emmm=Em-2;
Enmm=En-2;
Epm=Ep-1;
Epmm=Ep-2;
Euse=Edens{1}.Edensin;

%update the boundary for uhalfx and uhalfy

coef1=0.5*eust*dt/h./Euse;
uhalfz(Em,:,2:Epm)=uz(Em,:,2:Epm)+coef1(Em,:,2:Epm).*(p(Em,:,3:Ep)-p(Em,:,1:Epmm));
uhalfz(1,:,2:Epm)=uz(1,:,2:Epm)+coef1(1,:,2:Epm).*(p(1,:,3:Ep)-p(1,:,1:Epmm));
uhalfz(Em,:,Ep)=uz(Em,:,Ep)+coef1(Em,:,Ep).*(p(Em,:,2)-p(Em,:,Ep));
uhalfz(Em,:,1)=uz(Em,:,1)+coef1(Em,:,1).*(p(Em,:,2)-p(Em,:,Epm));
uhalfz(1,:,Ep)=uz(1,:,Ep)+coef1(1,:,Ep).*(p(1,:,2)-p(1,:,Epm));
uhalfz(1,:,1)=uz(1,:,1)+coef1(1,:,1).*(p(1,:,2)-p(1,:,Epm));

% uhalfy(1,:,:)=uy(1,:,:)+coef1(1,:,:).*(p(2,:,:)-p(1,:,:));
% uhalfy(Em,:,:)=uy(Em,:,:)+coef1(Em,:,:).*(p(Em,:,:)-p(Emm,:,:));
uhalfy=uy;

uhalfx(Em,1,:)=ux(Em,1,:)+coef1(Em,1,:).*(p(Em,2,:)-p(Em,Enm,:));
uhalfx(Em,En,:)=ux(Em,En,:)+coef1(Em,En,:).*(p(Em,2,:)-p(Em,Enm,:));
uhalfx(Em,2:Enm,:)=ux(Em,2:Enm,:)+coef1(Em,2:Enm,:).*(p(Em,3:En,:)-p(Em,1:Enmm,:));
uhalfx(1,2:Enm,:)=ux(1,2:Enm,:)+coef1(1,2:Enm,:).*(p(1,3:En,:)-p(1,1:Enmm,:));
uhalfx(1,1,:)=ux(1,1,:)+coef1(1,1,:).*(p(1,2,:)-p(1,Enm,:));
uhalfx(1,En,:)=ux(1,En,:)+coef1(1,En,:).*(p(1,2,:)-p(1,Enm,:));

uxux=uxm.^2;
uyuy=uym.^2;
uzuz=uzm.^2;
uzux=uzm.*uxm;
uyuz=uym.*uzm;
uyux=uxm.*uym;

rhsfixedx=zeros(Em,En,Ep);
rhsfixedy=rhsfixedx;
rhsfixedz=rhsfixedx;

rhsfixedx(2:Emm,2:Enm,2:Epm)=Euse(2:Emm,2:Enm,2:Epm).*(st*uxm(2:Emm,2:Enm,2:Epm)/dt-.25/h*(uxm(2:Emm,2:Enm,2:Epm).*...
            (uxm(2:Emm,3:En,2:Epm)-uxm(2:Emm,1:Enmm,2:Epm))+uzm(2:Emm,2:Enm,2:Epm).*(uxm(2:Emm,2:Enm,3:Ep)-...
            uxm(2:Emm,2:Enm,1:Epmm))+uym(2:Emm,2:Enm,2:Epm).*(uxm(3:Em,2:Enm,2:Epm)-uxm(1:Emmm,2:Enm,2:Epm))...
            +uxux(2:Emm,3:En,2:Epm)-uxux(2:Emm,1:Enmm,2:Epm)+uzux(2:Emm,2:Enm,3:Ep)-uzux(2:Emm,2:Enm,1:Epmm)...
            +uyux(3:Em,2:Enm,2:Epm)-uyux(1:Emmm,2:Enm,2:Epm)))+Efx(2:Emm,2:Enm,2:Epm);
rhsfixedy(2:Emm,2:Enm,2:Epm)=Euse(2:Emm,2:Enm,2:Epm).*(st*uym(2:Emm,2:Enm,2:Epm)/dt-.25/h*(uxm(2:Emm,2:Enm,2:Epm).*...
            (uym(2:Emm,3:En,2:Epm)-uym(2:Emm,1:Enmm,2:Epm))+uzm(2:Emm,2:Enm,2:Epm).*(uym(2:Emm,2:Enm,3:Ep)-...
            uym(2:Emm,2:Enm,1:Epmm))+uym(2:Emm,2:Enm,2:Epm).*(uym(3:Em,2:Enm,2:Epm)-uym(1:Emmm,2:Enm,2:Epm))...
            +uyux(2:Emm,3:En,2:Epm)-uyux(2:Emm,1:Enmm,2:Epm)+uyuz(2:Emm,2:Enm,3:Ep)-uyuz(2:Emm,2:Enm,1:Epmm)...
            +uyuy(3:Em,2:Enm,2:Epm)-uyuy(1:Emmm,2:Enm,2:Epm)))+Efy(2:Emm,2:Enm,2:Epm);
rhsfixedz(2:Emm,2:Enm,2:Epm)=Euse(2:Emm,2:Enm,2:Epm).*(st*uzm(2:Emm,2:Enm,2:Epm)/dt-.25/h*(uxm(2:Emm,2:Enm,2:Epm).*...
            (uzm(2:Emm,3:En,2:Epm)-uzm(2:Emm,1:Enmm,2:Epm))+uzm(2:Emm,2:Enm,2:Epm).*(uzm(2:Emm,2:Enm,3:Ep)-...
            uzm(2:Emm,2:Enm,1:Epmm))+uym(2:Emm,2:Enm,2:Epm).*(uzm(3:Em,2:Enm,2:Epm)-uzm(1:Emmm,2:Enm,2:Epm))...
            +uzux(2:Emm,3:En,2:Epm)-uzux(2:Emm,1:Enmm,2:Epm)+uzuz(2:Emm,2:Enm,3:Ep)-uzuz(2:Emm,2:Enm,1:Epmm)...
            +uyuz(3:Em,2:Enm,2:Epm)-uyuz(1:Emmm,2:Enm,2:Epm)))+Efz(2:Emm,2:Enm,2:Epm);
        
rhsfixedx(2:Emm,2:Enm,Ep)=Euse(2:Emm,2:Enm,Ep).*(st*uxm(2:Emm,2:Enm,Ep)/dt-0.25/h*(uxm(2:Emm,2:Enm,Ep).*...
    (uxm(2:Emm,3:En,Ep)-uxm(2:Emm,1:Enmm,Ep))+uzm(2:Emm,2:Enm,Ep).*(uxm(2:Emm,2:Enm,2)-...
    uxm(2:Emm,2:Enm,Epm))+uym(2:Emm,2:Enm,Ep).*(uxm(3:Em,2:Enm,Ep)-uxm(1:Emmm,2:Enm,Ep))...
    +uxux(2:Emm,3:En,Ep)-uxux(2:Emm,1:Enmm,Ep)+uzux(2:Emm,2:Enm,2)-uzux(2:Emm,2:Enm,Epm)...
    +uyux(3:Em,2:Enm,Ep)-uyux(1:Emmm,2:Enm,Ep)))+Efx(2:Emm,2:Enm,Ep);
rhsfixedx(2:Emm,2:Enm,1)=rhsfixedx(2:Emm,2:Enm,Ep);

rhsfixedy(2:Emm,2:Enm,Ep)=Euse(2:Emm,2:Enm,Ep).*(st*uym(2:Emm,2:Enm,Ep)/dt-.25/h*(uxm(2:Emm,2:Enm,Ep).*...
    (uym(2:Emm,3:En,Ep)-uym(2:Emm,1:Enmm,Ep))+uzm(2:Emm,2:Enm,Ep).*(uym(2:Emm,2:Enm,2)-...
    uym(2:Emm,2:Enm,Epm))+uym(2:Emm,2:Enm,Ep).*(uym(3:Em,2:Enm,Ep)-uym(1:Emmm,2:Enm,Ep))...
    +uyux(2:Emm,3:En,Ep)-uyux(2:Emm,1:Enmm,Ep)+uyuz(2:Emm,2:Enm,2)-uyuz(2:Emm,2:Enm,Epm)...
    +uyux(3:Em,2:Enm,Ep)-uyuy(1:Emmm,2:Enm,Ep)))+Efy(2:Emm,2:Enm,Ep);
rhsfixedy(2:Emm,2:Enm,1)=rhsfixedy(2:Emm,2:Enm,Ep);

rhsfixedz(2:Emm,2:Enm,Ep)=Euse(2:Emm,2:Enm,Ep).*(st*uzm(2:Emm,2:Enm,Ep)/dt-.25/h*(uxm(2:Emm,2:Enm,Ep).*...
    (uzm(2:Emm,3:En,Ep)-uzm(2:Emm,1:Enmm,Ep))+uzm(2:Emm,2:Enm,Ep).*(uzm(2:Emm,2:Enm,2)-...
    uzm(2:Emm,2:Enm,Epm))+uym(2:Emm,2:Enm,Ep).*(uzm(3:Em,2:Enm, Ep)-uzm(1:Emmm,2:Enm,Ep))...
    +uzux(2:Emm,3:En,Ep)-uzux(2:Emm,1:Enmm,Ep)+uzuz(2:Emm,2:Enm,2)-uzuz(2:Emm,2:Enm,Epm)...
    +uyuz(3:Em,2:Enm,Ep)-uyuz(1:Emmm,2:Enm,Ep)))+Efz(2:Emm,2:Enm,Ep);
rhsfixedz(2:Emm,2:Enm,1)=rhsfixedz(2:Emm,2:Enm,Ep);

rhsfixedx(2:Emm,En,2:Epm)=Euse(2:Emm,En,2:Epm).*(st*uxm(2:Emm,En,2:Epm)/dt-.25/h*(uxm(2:Emm,En,2:Epm).*...
            (uxm(2:Emm,2,2:Epm)-uxm(2:Emm,Enm,2:Epm))+uzm(2:Emm,En,2:Epm).*(uxm(2:Emm,En,3:Ep)-...
            uxm(2:Emm,En,1:Epmm))+uym(2:Emm,En,2:Epm).*(uxm(3:Em,En,2:Epm)-uxm(1:Emmm,En,2:Epm))...
            +uxux(2:Emm,2,2:Epm)-uxux(2:Emm,Enm,2:Epm)+uzux(2:Emm,En,3:Ep)-uzux(2:Emm,En,1:Epmm)...
            +uyux(3:Em,En,2:Epm)-uyux(1:Emmm,En,2:Epm)))+Efx(2:Emm,En,2:Epm);
rhsfixedx(2:Emm,1,2:Epm)=rhsfixedx(2:Emm,En,2:Epm);

rhsfixedy(2:Emm,En,2:Epm)=Euse(2:Emm,En,2:Epm).*(st*uym(2:Emm,En,2:Epm)/dt-.25/h*(uxm(2:Emm,En,2:Epm).*...
            (uym(2:Emm,2,2:Epm)-uym(2:Emm,Enm,2:Epm))+uzm(2:Emm,En,2:Epm).*(uym(2:Emm,En,3:Ep)-...
            uym(2:Emm,En,1:Epmm))+uym(2:Emm,En,2:Epm).*(uym(3:Em,En,2:Epm)-uym(1:Emmm,En,2:Epm))...
            +uyux(2:Emm,2,2:Epm)-uyux(2:Emm,Enm,2:Epm)+uyuz(2:Emm,En,3:Ep)-uyuz(2:Emm,En,1:Epmm)...
            +uyuy(3:Em,En,2:Epm)-uyuy(1:Emmm,En,2:Epm)))+Efy(2:Emm,En,2:Epm);
rhsfixedy(2:Emm,1,2:Epm)=rhsfixedy(2:Emm,En,2:Epm);

rhsfixedz(2:Emm,En,2:Epm)=Euse(2:Emm,En,2:Epm).*(st*uzm(2:Emm,En,2:Epm)/dt-.25/h*(uxm(2:Emm,En,2:Epm).*...
            (uzm(2:Emm,2,2:Epm)-uzm(2:Emm,Enm,2:Epm))+uzm(2:Emm,En,2:Epm).*(uzm(2:Emm,En,3:Ep)-...
            uzm(2:Emm,En,1:Epmm))+uym(2:Emm,En,2:Epm).*(uzm(3:Em,En,2:Epm)-uzm(1:Emmm,En,2:Epm))...
            +uzux(2:Emm,2,2:Epm)-uzux(2:Emm,Enm,2:Epm)+uzuz(2:Emm,En,3:Ep)-uzuz(2:Emm,En,1:Epmm)...
            +uyuz(3:Em,En,2:Epm)-uyuz(1:Emmm,En,2:Epm)))+Efz(2:Emm,En,2:Epm);
rhsfixedz(2:Emm,1,2:Epm)=rhsfixedz(2:Emm,En,2:Epm);

rhsfixedx(2:Emm,En,Ep)=Euse(2:Emm,En,Ep).*(st*uxm(2:Emm,En,Ep)/dt-.25/h*(uxm(2:Emm,En,Ep).*...
            (uxm(2:Emm,2,Ep)-uxm(2:Emm,Enm,Ep))+uzm(2:Emm,En,Ep).*(uxm(2:Emm,En,2)-...
            uxm(2:Emm,En,Epm))+uym(2:Emm,En,Ep).*(uxm(3:Em,En,Ep)-uxm(1:Emmm,En,Ep))...
            +uxux(2:Emm,2,Ep)-uxux(2:Emm,Enm,Ep)+uzux(2:Emm,En,2)-uzux(2:Emm,En,Epm)...
            +uyux(3:Em,En,Ep)-uyux(1:Emmm,En,Ep)))+Efx(2:Emm,En,Ep);
rhsfixedx(2:Emm,1,Ep)=rhsfixedx(2:Emm,En,Ep);
rhsfixedx(2:Emm,1,1)=rhsfixedx(2:Emm,En,Ep);
rhsfixedx(2:Emm,En,1)=rhsfixedx(2:Emm,En,Ep);

rhsfixedy(2:Emm,En,Ep)=Euse(2:Emm,En,Ep).*(st*uym(2:Emm,En,Ep)/dt-.25/h*(uxm(2:Emm,En,Ep).*...
            (uym(2:Emm,2,Ep)-uym(2:Emm,Enm,Ep))+uzm(2:Emm,En,Ep).*(uym(2:Emm,En,2)-...
            uym(2:Emm,En,Epm))+uym(2:Emm,En,Ep).*(uym(3:Em,En,Ep)-uym(1:Emmm,En,Ep))...
            +uyux(2:Emm,2,Ep)-uyux(2:Emm,Enm,Ep)+uyuz(2:Emm,En,2)-uyuz(2:Emm,En,Epm)...
            +uyuy(3:Em,En,Ep)-uyuy(1:Emmm,En,Ep)))+Efy(2:Emm,En,Ep);
rhsfixedy(2:Emm,1,Ep)=rhsfixedy(2:Emm,En,Ep);
rhsfixedy(2:Emm,1,1)=rhsfixedy(2:Emm,En,Ep);
rhsfixedy(2:Emm,En,1)=rhsfixedy(2:Emm,En,Ep);

rhsfixedz(2:Emm,En,Ep)=Euse(2:Emm,En,Ep).*(st*uzm(2:Emm,En,Ep)/dt-.25/h*(uxm(2:Emm,En,Ep).*...
            (uzm(2:Emm,2,Ep)-uzm(2:Emm,Enm,Ep))+uzm(2:Emm,En,Ep).*(uzm(2:Emm,En,2)-...
            uzm(2:Emm,En,Epm))+uym(2:Emm,En,Ep).*(uzm(3:Em,En,Ep)-uzm(1:Emmm,En,Ep))...
            +uzux(2:Emm,2,Ep)-uzux(2:Emm,Enm,Ep)+uzuz(2:Emm,En,2)-uzuz(2:Emm,En,Epm)...
            +uyuz(3:Em,En,Ep)-uyuz(1:Emmm,En,Ep)))+Efz(2:Emm,En,Ep);
rhsfixedz(2:Emm,1,Ep)=rhsfixedz(2:Emm,En,Ep);
rhsfixedz(2:Emm,1,1)=rhsfixedz(2:Emm,En,Ep);
rhsfixedz(2:Emm,En,1)=rhsfixedz(2:Emm,En,Ep);

rhsfixedx(:,:,1)=rhsfixedx(:,:,Ep); rhsfixedx(:,1,:)=rhsfixedx(:,En,:);
rhsfixedy(:,:,1)=rhsfixedy(:,:,Ep); rhsfixedy(:,1,:)=rhsfixedy(:,En,:); 
rhsfixedz(:,:,1)=rhsfixedz(:,:,Ep); rhsfixedz(:,1,:)=rhsfixedz(:,En,:);


resx=zeros(Em,En,Ep);
resy=zeros(Em,En,Ep);
resz=zeros(Em,En,Ep);
[vxsave,vysave,vzsave]=fmvVEL3Dper2(uhalfx,uhalfy,uhalfz,Re*rhsfixedx,Re*rhsfixedy,Re*rhsfixedz,h,Edens,levels,vcoef);
[resx,resy,resz]=residualvel3Dper2(vxsave,vysave,vzsave,Re*rhsfixedx,Re*rhsfixedy,Re*rhsfixedz,vcoef{1}, resx, resy, resz);

vxerr=max(max(max(abs(resx))));
vyerr=max(max(max(abs(resy))));
vzerr=max(max(max(abs(resz))));

[zx,zy,zz]=fmvVEL3Dper2(zeros(Em,En,Ep),zeros(Em,En,Ep),zeros(Em,En,Ep),resx,resy,resz,h,Edens,levels,vcoef);

p0x=zx;
p0y=zy;
p0z=zz;


vxnew=vxsave;
vxerrsave=vxerr;

vynew=vysave;
vyerrsave=vyerr;

vznew=vzsave;
vzerrsave=vzerr;

rho=sum(sum(sum(resx(2:Emm,2:En,1:Epm).*zx(2:Emm,2:En,1:Epm)+resy(2:Emm,2:En,1:Epm).*zy(2:Emm,2:En,1:Epm)...
    +resz(2:Emm,2:En,1:Epm).*zz(2:Emm,2:En,1:Epm))));

count1=0;
c=0;

while (vxerr>allowerr || vyerr>allowerr || vzerr>allowerr) && count1<2 && c<20
    
    [Ap0x,Ap0y,Ap0z]=AonVecVEL3Dper2(p0x,p0y,p0z,vcoef{1}); %h,dt,Edens{1},viscmatmid{1},viscmat{1},st,Re);
    alp=rho/sum(sum(sum(p0x(2:Emm,2:En,1:Epm).*Ap0x(2:Emm,2:En,1:Epm)+p0y(2:Emm,2:En,1:Epm).*Ap0y(2:Emm,2:En,1:Epm)...
        +p0z(2:Emm,2:En,1:Epm).*Ap0z(2:Emm,2:En,1:Epm))));
    vxnew=vxnew+alp*p0x;
    vynew=vynew+alp*p0y;
    vznew=vznew+alp*p0z;
    
  
    [resx,resy,resz]=residualvel3Dper2(vxnew,vynew,vznew,Re*rhsfixedx,Re*rhsfixedy,Re*rhsfixedz,vcoef{1}, resx,resy,resz);
    vxerr=max(max(max(abs(resx))));
    vyerr=max(max(max(abs(resy))));
    vzerr=max(max(max(abs(resz))));
    
    [zx,zy,zz]=mvVEL3Dper2(zeros(Em,En,Ep),zeros(Em,En,Ep),zeros(Em,En,Ep),resx,resy,resz,h,Edens,levels,vcoef);

    for i=1:3
        [zx,zy,zz]=mvVEL3Dper2(zx,zy,zz,resx,resy,resz,h,Edens,levels,vcoef);
    end
      rhonew=sum(sum(sum(resx(2:Emm,2:En,1:Epm).*zx(2:Emm,2:En,1:Epm)+resy(2:Emm,2:En,1:Epm).*zy(2:Emm,2:En,1:Epm)...
          +resz(2:Emm,2:En,1:Epm).*zz(2:Emm,2:En,1:Epm))));
    bet=rhonew/rho;
    p0x=zx+bet*p0x;
    p0y=zy+bet*p0y;
    p0z=zz+bet*p0z;
    rho=rhonew;
    
    
  %     max(max(abs(rhsfixedx)))
%     vxerrscl=vxerr/max(max(abs(rhsfixedx)))*max(max(abs(vxnew)))
%     
    vxdiv=vxerrsave/vxerr;
    vydiv=vyerrsave/vyerr;
    vzdiv=vzerrsave/vzerr;
    
     if vxerr<vxerrsave || vyerr<vyerrsave || vzerr<vzerrsave
        
        if vxerr<vxerrsave
        vxsave=vxnew;
        vxerrsave=vxerr;
        end
        
        if vyerr<vyerrsave        
        vysave=vynew;
        vyerrsave=vyerr;
        end
        
        if vzerr<vzerrsave        
        vzsave=vznew;
        vzerrsave=vzerr;
        end
    else
    count1=count1+1;
    end
    
    c=c+1;
end

if count1>=2 && vzerrsave>10^-2
    vxsave=0;
    vysave=0;
    vzsave=0;
end
