function [ux,uy,uz]=DS_Vel2(ux,uy,uz,p,eust,dt,h,Edens,Em,En,Ep,uhalfx,uhalfy,uhalfz)
    Emm=Em-1; Enm=En-1; Epm=Ep-1; Epmm=Ep-2; Emmm=Em-2; Enmm=En-2;
    Edens1=Edens{1}.Edensin;

    edenscoefx=eust*dt/(h)*(1./(Edens1(2:Emm,1:Enm,2:Ep)+Edens1(2:Emm,2:En,2:Ep)));
    edenscoefy=eust*dt/(h)*(1./(Edens1(2:Emm,2:En,2:Ep)+Edens1(3:Em,2:En,2:Ep)));
    edenscoefz=eust*dt/(h)*(1./(Edens1(2:Emm,2:En,1:Epm)+Edens1(2:Emm,2:En,2:Ep)));
    edenscoef=eust*dt/h*0.5./Edens1;
    
    ux(2:Emm,2:Enm,1:Ep)=(p(2:Emm,1:Enmm,1:Ep)-...
        p(2:Emm,3:En,1:Ep)).*edenscoef(2:Emm,2:Enm,1:Ep)+uhalfx(2:Emm,2:Enm,1:Ep);
   
    ux(2:Emm,1,1:Ep)=(p(2:Emm,Enm,1:Ep)-...
        p(2:Emm,2,1:Ep)).*edenscoef(2:Emm,1,1:Ep)+uhalfx(2:Emm,1,1:Ep);
    
    ux(2:Emm,En,1:Ep)=ux(2:Emm,1,1:Ep);
    
    uy(2:Emm,1:En,1:Ep)=(p(1:Emmm,1:En,1:Ep)-...
        p(3:Em,1:En,1:Ep)).*edenscoef(2:Emm,:,:)+uhalfy(2:Emm,1:En,1:Ep);
    
    uz(2:Emm,1:En,2:Epm)=(p(2:Emm,1:En,1:Epmm)-...
        p(2:Emm,1:En,3:Ep)).*edenscoef(2:Emm,:,2:Epm)+uhalfz(2:Emm,1:En,2:Epm);
    uz(2:Emm,1:En,Ep)=(p(2:Emm,1:En,Epm)-...
        p(2:Emm,1:En,2)).*edenscoef(2:Emm,:,Ep)+uhalfz(2:Emm,1:En,Ep);
    uz(2:Emm,1:En,1)=uz(2:Emm,1:En,Ep);