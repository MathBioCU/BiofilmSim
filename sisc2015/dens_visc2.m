function [Edens,viscmat, viscmatmid]=dens_visc2(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,ylength, xlength,visc,addlvisc,levelsP,Em,En,Ep,r)
Edens{1}.dx=h;
Edens{1}.Edensin=transferLtoEdens3Dper_e3a(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength,r);
viscmat{1}=transferLtoEvisc3Dper_e3a(h,x,y,z,X,visc,addlvisc,numOfnonzero,1,zlength,xlength,r);
        Edens{1}.Edensmidlr=(Edens{1}.Edensin(:,2:En,:)+Edens{1}.Edensin(:,1:En-1,:))/2;
        Edens{1}.Edensmidud=(Edens{1}.Edensin(2:Em,:,:)+Edens{1}.Edensin(1:Em-1,:,:))/2;
        Edens{1}.Edensmidfb=(Edens{1}.Edensin(:,:,2:Ep)+Edens{1}.Edensin(:,:,1:Ep-1))/2;
        
        Edens{1}.iEdensin=1./Edens{1}.Edensin;
        Edens{1}.iEdensmidlr=(Edens{1}.iEdensin(:,2:En,:)+Edens{1}.iEdensin(:,1:En-1,:))/2;
        Edens{1}.iEdensmidud=(Edens{1}.iEdensin(2:Em,:,:)+Edens{1}.iEdensin(1:Em-1,:,:))/2;
        Edens{1}.iEdensmidfb=(Edens{1}.iEdensin(:,:,2:Ep)+Edens{1}.iEdensin(:,:,1:Ep-1))/2;        
       
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

for c5=2:levelsP
    clear x2h y2h z2h dx hytemp yvectemp twotopow Em1 En1 Ep1 
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
           
        viscmat{c5}=restricthto2h3DVper2(viscmat{c5-1});
        viscmatmid{c5}.lr=(viscmat{c5}(:,2:En1,:)+viscmat{c5}(:,1:En1-1,:))/2;
        viscmatmid{c5}.ud=(viscmat{c5}(2:Em1,:,:)+viscmat{c5}(1:Em1-1,:,:))/2;
        viscmatmid{c5}.fb=(viscmat{c5}(:,:,2:Ep1)+viscmat{c5}(:,:,1:Ep1-1))/2;
        
        Edens{c5}.Edensin=restricthto2h3Dper2(Edens{c5-1}.Edensin);
        Edens{c5}.Edensmidlr=(Edens{c5}.Edensin(:,2:En1,:)+Edens{c5}.Edensin(:,1:En1-1,:))/2;
        Edens{c5}.Edensmidud=(Edens{c5}.Edensin(2:Em1,:,:)+Edens{c5}.Edensin(1:Em1-1,:,:))/2;
        Edens{c5}.Edensmidfb=(Edens{c5}.Edensin(:,:,2:Ep1)+Edens{c5}.Edensin(:,:,1:Ep1-1))/2;
        
        
        Edens{c5}.iEdensin=1./Edens{c5}.Edensin;
        Edens{c5}.iEdensmidlr=(Edens{c5}.iEdensin(:,2:En1,:)+Edens{c5}.iEdensin(:,1:En1-1,:))/2;
        Edens{c5}.iEdensmidud=(Edens{c5}.iEdensin(2:Em1,:,:)+Edens{c5}.iEdensin(1:Em1-1,:,:))/2;
        Edens{c5}.iEdensmidfb=(Edens{c5}.iEdensin(:,:,2:Ep1)+Edens{c5}.iEdensin(:,:,1:Ep1-1))/2;      
       
    
end