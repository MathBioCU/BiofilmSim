function [ux,uy,uz]=DS_Vel(ux,uy,uz,p,eust,dt,h,Edens,Em,En,Ep,uhalfx,uhalfy,uhalfz)
    Emm=Em-1; Enm=En-1; Epm=Ep-1; Epmm=Ep-2; Emmm=Em-2; Enmm=En-2;
    Edens1=Edens{1}.Edensin;

    edenscoefx=2*eust*dt/(h)*(1./(Edens1(2:Emm,1:Enm,2:Ep)+Edens1(2:Emm,2:En,2:Ep)));
    edenscoefy=2*eust*dt/(h)*(1./(Edens1(2:Emm,2:En,2:Ep)+Edens1(3:Em,2:En,2:Ep)));
    edenscoefz=2*eust*dt/(h)*(1./(Edens1(2:Emm,2:En,1:Epm)+Edens1(2:Emm,2:En,2:Ep)));
  
    ux(2:Emm,1:Enm,2:Ep)=(p(2:Emm,1:Enm,2:Ep)-...
        p(2:Emm,2:En,2:Ep)).*edenscoefx(:,1:Enm,1:Epm)+uhalfx(2:Emm,1:Enm,2:Ep);
    
    ux(2:Emm,1:Enm,Ep)=(p(2:Emm,1:Enm,Ep)-...
        p(2:Emm,2:En,Ep)).*edenscoefx(:,1:Enm,Epm)+uhalfx(2:Emm,1:Enm,Ep);
    ux(2:Emm,1:Enm,1)=ux(2:Emm,1:Enm,Ep);
   
    ux(2:Emm,En,2:Epm)=ux(2:Emm,1,2:Epm);
    
%     ux(2:Emm,1,Ep)=(p(2:Emm,Enm,Ep)-p(2:Emm,En,Ep)).*edenscoefx(:,Enm,Epm)+uhalfx(2:Emm,En,Ep);
    ux(2:Emm,1,1)=ux(2:Emm,1,Ep);
    ux(2:Emm,En,Ep)=ux(2:Emm,1,Ep);
    ux(2:Emm,En,1)=ux(2:Emm,1,Ep);
    
    
    uy(2:Emm,2:En,2:Epm)=(p(2:Emm,2:En,2:Epm)-...
        p(3:Em,2:En,2:Epm)).*edenscoefy(:,:,1:Epmm)+uhalfy(2:Emm,2:En,2:Epm);
    uy(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Ep)-...
        p(3:Em,2:En,Ep)).*edenscoefy(:,:,Epm)+uhalfy(2:Emm,2:En,Ep);
    uy(2:Emm,2:En,1)=uy(2:Emm,2:En,Ep);
    uy(2:Emm,1,:)=uy(2:Emm,En,:);
    
    uz(2:Emm,2:En,1:Epm)=(p(2:Emm,2:En,1:Epm)-...
        p(2:Emm,2:En,2:Ep)).*edenscoefz(:,:,1:Epm)+uhalfz(2:Emm,2:En,1:Epm);
    uz(2:Emm,2:En,Ep)=(p(2:Emm,2:En,Ep)-...
        p(2:Emm,2:En,2)).*edenscoefz(:,:,Epm)+uhalfz(2:Emm,2:En,Ep);
    uz(2:Emm,2:En,Ep)=uz(2:Emm,2:En,1);
    uz(2:Emm,1,:)=uz(2:Emm,En,:);