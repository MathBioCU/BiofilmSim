function [vcoef,pcoef]=compute_operator(Edens,levelsP,Re,st,viscmat, viscmatmid,dt)
for c13=1:levelsP
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
       
    pcoef{c13}.coef=(1*Edens{c13}.iEdensmidlr(2:tEmm,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidlr(2:tEmm,1:tEnmm,2:tEpm)+...
        1*Edens{c13}.iEdensmidud(2:tEmm,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidud(1:tEmmm,2:tEnm,2:tEpm)+...
        1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    
    %Neumann boundary not at periodic ends
    pcoef{c13}.coefboundB=(1*Edens{c13}.iEdensmidlr(1,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidlr(1,1:tEnmm,2:tEpm)+...
        2*Edens{c13}.iEdensmidud(1,2:tEnm,2:tEpm)+...
        1*Edens{c13}.iEdensmidfb(1,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidfb(1,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT=(1*Edens{c13}.iEdensmidlr(tEm,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidlr(tEm,1:tEnmm,2:tEpm)+...
        2*Edens{c13}.iEdensmidud(tEmm,2:tEnm,2:tEpm)+...
        1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,2:tEpm)+1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,1:tEpmm))/Edens{c13}.dx^2;
    
    %formula for non-edge points on periodic ends
    pcoef{c13}.coef(1:tEmmm,1:tEnmm,tEpm)=(1*Edens{c13}.iEdensmidlr(2:tEmm,2:tEnm,1)+1*Edens{c13}.iEdensmidlr(2:tEmm,1:tEnmm,1)+... %adjust indices to get right size
        1*Edens{c13}.iEdensmidud(2:tEmm,2:tEnm,1)+1*Edens{c13}.iEdensmidud(1:tEmmm,2:tEnm,1)+...
        1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,1)+1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coef(1:tEmmm,tEnm,1:tEpmm)=(1*Edens{c13}.iEdensmidlr(2:tEmm,1,2:tEpm)+1*Edens{c13}.iEdensmidlr(2:tEmm,tEnm,2:tEpm)+...
        1*Edens{c13}.iEdensmidud(2:tEmm,1,2:tEpm)+1*Edens{c13}.iEdensmidud(1:tEmmm,1,2:tEpm)+...
        1*Edens{c13}.iEdensmidfb(2:tEmm,1,2:tEpm)+1*Edens{c13}.iEdensmidfb(2:tEmm,1,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coef(1:tEmmm,tEnm,tEpm)=(1*Edens{c13}.iEdensmidlr(2:tEmm,1,1)+1*Edens{c13}.iEdensmidlr(2:tEmm,tEnm,1)+...
        1*Edens{c13}.iEdensmidud(2:tEmm,1,1)+1*Edens{c13}.iEdensmidud(1:tEmmm,1,1)+...
        1*Edens{c13}.iEdensmidfb(2:tEmm,1,1)+1*Edens{c13}.iEdensmidfb(2:tEmm,1,tEpm))/Edens{c13}.dx^2;
    
    %edge-point on periodic boundary
    pcoef{c13}.coefboundB(1,1:tEnmm,tEpm)=(1*Edens{c13}.iEdensmidlr(1,2:tEnm,tEp)+1*Edens{c13}.iEdensmidlr(1,1:tEnmm,tEp)+...
        2*Edens{c13}.iEdensmidud(1,2:tEnm,tEp)+...
        1*Edens{c13}.iEdensmidfb(1,2:tEnm,1)+1*Edens{c13}.iEdensmidfb(1,2:tEnm,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundB(1,tEnm,1:tEpmm)=(1*Edens{c13}.iEdensmidlr(1,1,2:tEpm)+1*Edens{c13}.iEdensmidlr(1,tEnm,2:tEpm)+...
        2*Edens{c13}.iEdensmidud(1,tEn,2:tEpm)+...
        1*Edens{c13}.iEdensmidfb(1,tEn,2:tEpm)+1*Edens{c13}.iEdensmidfb(1,tEn,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT(1,1:tEnmm,tEpm)=(1*Edens{c13}.iEdensmidlr(tEm,2:tEnm,tEp)+1*Edens{c13}.iEdensmidlr(tEm,1:tEnmm,tEp)+...
        2*Edens{c13}.iEdensmidud(tEmm,2:tEnm,tEp)+...
        1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,1)+1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT(1,tEnm,1:tEpmm)=(1*Edens{c13}.iEdensmidlr(tEm,1,2:tEpm)+1*Edens{c13}.iEdensmidlr(tEm,tEnm,2:tEpm)+...
        2*Edens{c13}.iEdensmidud(tEmm,tEn,2:tEpm)+...
        1*Edens{c13}.iEdensmidfb(tEm,tEn,2:tEpm)+1*Edens{c13}.iEdensmidfb(tEm,tEn,1:tEpmm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundB(1,tEnm,tEpm)=(1*Edens{c13}.iEdensmidlr(1,1,tEp)+1*Edens{c13}.iEdensmidlr(1,tEnm,tEp)+...
        2*Edens{c13}.iEdensmidud(1,tEn,tEp)+...
        1*Edens{c13}.iEdensmidfb(1,tEn,1)+1*Edens{c13}.iEdensmidfb(1,tEn,tEpm))/Edens{c13}.dx^2;
    pcoef{c13}.coefboundT(1,tEnm,tEpm)=(1*Edens{c13}.iEdensmidlr(tEm,1,tEp)+1*Edens{c13}.iEdensmidlr(tEm,tEnm,tEp)+...
        2*Edens{c13}.iEdensmidud(tEmm,tEn,tEp)+...
        1*Edens{c13}.iEdensmidfb(tEm,tEn,1)+1*Edens{c13}.iEdensmidfb(tEm,tEn,tEpm))/Edens{c13}.dx^2;

       
    pcoef{c13}.rescoefpp1=1*Edens{c13}.iEdensmidlr(2:tEmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp1(1:tEmmm,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidlr(2:tEmm,1,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp1boundB=1*Edens{c13}.iEdensmidlr(1,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp1boundB(1,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidlr(1,1,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp1boundT=1*Edens{c13}.iEdensmidlr(tEm,2:tEnm,2:tEp)/Edens{c13}.dx^2; %
    pcoef{c13}.rescoefpp1boundT(1,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidlr(tEm,1,2:tEp)/Edens{c13}.dx^2;
 
    pcoef{c13}.rescoefpm1=1*Edens{c13}.iEdensmidlr(2:tEmm,1:tEnmm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1(1:tEmmm,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidlr(2:tEmm,tEnm,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpm1boundB=1*Edens{c13}.iEdensmidlr(1,1:tEnmm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1boundB(1,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidlr(1,tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1boundT=1*Edens{c13}.iEdensmidlr(tEm,1:tEnmm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm1boundT(1,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidlr(tEm,tEnm,2:tEp)/Edens{c13}.dx^2;


    pcoef{c13}.rescoefpp2=1*Edens{c13}.iEdensmidud(2:tEmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp2(1:tEmmm,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidud(2:tEmm,tEn,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp2boundB=2*Edens{c13}.iEdensmidud(1,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp2boundB(1,tEnm,1:tEpm)=2*Edens{c13}.iEdensmidud(1,tEn,2:tEp)/Edens{c13}.dx^2;

    pcoef{c13}.rescoefpm2=1*Edens{c13}.iEdensmidud(1:tEmmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm2(1:tEmmm,tEnm,1:tEpm)=1*Edens{c13}.iEdensmidud(1:tEmmm,tEn,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpm2boundT=2*Edens{c13}.iEdensmidud(tEmm,2:tEnm,2:tEp)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm2boundT(1,tEnm,1:tEpm)=2*Edens{c13}.iEdensmidud(tEmm,tEn,2:tEp)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp3=1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3(1:tEmmm,1:tEnmm,tEpm)=1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,1)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3(1:tEmmm,tEnm,1:tEpmm)=1*Edens{c13}.iEdensmidfb(2:tEmm,tEn,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3(1:tEmmm,tEnm,tEpm)=1*Edens{c13}.iEdensmidfb(2:tEmm,tEn,1)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp3boundB=1*Edens{c13}.iEdensmidfb(1,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundB(1,1:tEnmm,tEpm)=1*Edens{c13}.iEdensmidfb(1,2:tEnm,1)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundB(1,tEnm,1:tEpmm)=1*Edens{c13}.iEdensmidfb(1,tEn,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundB(1,tEnm,tEpm)=1*Edens{c13}.iEdensmidfb(1,tEn,1)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpp3boundT=1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundT(1,1:tEnmm,tEpm)=1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,1)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundT(1,tEnm,1:tEpmm)=1*Edens{c13}.iEdensmidfb(tEm,tEn,2:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpp3boundT(1,tEnm,tEpm)=1*Edens{c13}.iEdensmidfb(tEm,tEn,1)/Edens{c13}.dx^2;

    pcoef{c13}.rescoefpm3=1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,1:tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3(1:tEmmm,1:tEnmm,tEpm)=1*Edens{c13}.iEdensmidfb(2:tEmm,2:tEnm,tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3(1:tEmmm,tEnm,1:tEpmm)=1*Edens{c13}.iEdensmidfb(2:tEmm,tEn,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3(1:tEmmm,tEnm,tEpm)=1*Edens{c13}.iEdensmidfb(2:tEmm,tEn,tEpm)/Edens{c13}.dx^2;
    
    pcoef{c13}.rescoefpm3boundB=1*Edens{c13}.iEdensmidfb(1,2:tEnm,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundB(1,:,tEpm)=1*Edens{c13}.iEdensmidfb(1,2:tEnm,tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundB(1,tEnm,1:tEpmm)=1*Edens{c13}.iEdensmidfb(1,tEn,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundB(1,tEnm,tEpm)=1*Edens{c13}.iEdensmidfb(1,tEn,tEpm)/Edens{c13}.dx^2;

    pcoef{c13}.rescoefpm3boundT=1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundT(1,:,tEpm)=1*Edens{c13}.iEdensmidfb(tEm,2:tEnm,tEpm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundT(1,tEnm,1:tEpmm)=1*Edens{c13}.iEdensmidfb(tEm,tEn,1:tEpmm)/Edens{c13}.dx^2;
    pcoef{c13}.rescoefpm3boundT(1,tEnm,tEpm)=1*Edens{c13}.iEdensmidfb(tEm,tEn,tEpm)/Edens{c13}.dx^2;
    
    
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
