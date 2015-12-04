function [Edens,viscmat,viscmatmid]=DensVisc(addlvisc,addldens,h,x,y,z,X,visc,numOfnonzero,xlength,zlength,Em,En,Ep,d0mean,initdensity,levelsV,viscmat,Edens)
if addlvisc>0 || addldens>0 || h>0
        viscmat{1}=transferLtoEvisc3Dper_e4a(h,x,y,z,X,visc,addlvisc,numOfnonzero,1,zlength,xlength);
        viscmatmid{1}.lr=(viscmat{1}(:,2:En,:)+viscmat{1}(:,1:En-1,:))/2;
        viscmatmid{1}.ud=(viscmat{1}(2:Em,:,:)+viscmat{1}(1:Em-1,:,:))/2;
        viscmatmid{1}.fb=(viscmat{1}(:,:,2:Ep)+viscmat{1}(:,:,1:Ep-1))/2;
        
     Edens{1}.Edensin=transferLtoEdens3Dper_e2a(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength);
     Edens{1}.Edensmidlr=(Edens{1}.Edensin(:,2:En,:)+Edens{1}.Edensin(:,1:En-1,:))/2;
     Edens{1}.Edensmidud=(Edens{1}.Edensin(2:Em,:,:)+Edens{1}.Edensin(1:Em-1,:,:))/2;
     Edens{1}.Edensmidfb=(Edens{1}.Edensin(:,:,2:Ep)+Edens{1}.Edensin(:,:,1:Ep-1))/2;

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