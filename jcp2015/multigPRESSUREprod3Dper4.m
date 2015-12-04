function [psave,perrsave]=multigPRESSUREprod3Dper4(dt,h,p,uhalfx,uhalfy,uhalfz,Edens,Em,En,Ep,levels,pcoef,eust)

Emm=Em-1;
Enm=En-1;
Emmm=Em-2;
Enmm=En-2;
Epm=Ep-1;
Epmm=Ep-2;

% relaxtimes=35;
% allowerr=10^-20;
allowerr=10^-9;
% allowerr=10^-8;
% En=En+1;

Euse=Edens{1}.Edensin;
rhsfixed=zeros(Em,En,Ep);

% uhalfx=[uhalfx,uhalfx(:,2)];
% uhalfy=[uhalfy,uhalfy(:,2)];

rhsfixed(2:Emm,2:Enm,2:Epm)=(uhalfx(2:Emm,2:Enm,2:Epm)-uhalfx(2:Emm,1:Enmm,2:Epm)+uhalfy(2:Emm,2:Enm,2:Epm)-...
    uhalfy(1:Emmm,2:Enm,2:Epm)+uhalfz(2:Emm,2:Enm,2:Epm)-uhalfz(2:Emm,2:Enm,1:Epmm))/dt/eust;

rhsfixed(2:Emm,2:Enm,1)=(uhalfx(2:Emm,2:Enm,1)-uhalfx(2:Emm,1:Enmm,1)+uhalfy(2:Emm,2:Enm,1)-...
    uhalfy(1:Emmm,2:Enm,1)+uhalfz(2:Emm,2:Enm,Ep)-uhalfz(2:Emm,2:Enm,Epm))/dt/eust;  %periodic condition
rhsfixed(2:Emm,2:Enm,Ep)=rhsfixed(2:Emm,2:Enm,1);

rhsfixed(1,2:Enm,2:Epm)=(uhalfx(1,2:Enm,2:Epm)-uhalfx(1,1:Enmm,2:Epm)+...
    uhalfz(1,2:Enm,2:Epm)-uhalfz(1,2:Enm,1:Epmm))/dt/eust; %bottom boundary

rhsfixed(1,2:Enm,Ep)=(uhalfx(1,2:Enm,Ep)-uhalfx(1,1:Enmm,Ep)...
    +uhalfz(1,2:Enm,1)-uhalfz(1,2:Enm,Epm))/dt/eust; %bottom boundary

rhsfixed(Em,2:Enm,2:Epm)=(uhalfx(Em,2:Enm,2:Epm)-uhalfx(Em,1:Enmm,2:Epm)...
    +uhalfz(Em,2:Enm,2:Epm)-uhalfz(Em,2:Enm,1:Epmm))/dt/eust; %top boundary

rhsfixed(Em,2:Enm,Ep)=(uhalfx(Em,2:Enm,Ep)-uhalfx(Em,1:Enmm,Ep)...
    +uhalfz(Em,2:Enm,1)-uhalfz(Em,2:Enm,Epm))/dt/eust; %top boundary

rhsfixed(2:Emm,1,2:Epm)=(uhalfx(2:Emm,1,2:Epm)-uhalfx(2:Emm,Enm,2:Epm)+uhalfy(2:Emm,En,2:Epm)-...
    uhalfy(1:Emmm,En,2:Epm)+uhalfz(2:Emm,En,2:Epm)-uhalfz(2:Emm,En,1:Epmm))/dt/eust; %LR periodic

rhsfixed(2:Emm,1,Ep)=(uhalfx(2:Emm,1,Ep)-uhalfx(2:Emm,Enm,Ep)+uhalfy(2:Emm,En,Ep)-...
    uhalfy(1:Emmm,En,Ep)+uhalfz(2:Emm,En,1)-uhalfz(2:Emm,En,Epm))/dt/eust; %left boundary

rhsfixed(2:Emm,En,Ep)=rhsfixed(2:Emm,1,Ep);
rhsfixed(2:Emm,1,1)=rhsfixed(2:Emm,1,Ep);
rhsfixed(2:Emm,En,1)=rhsfixed(2:Emm,1,Ep);

rhsfixed(1,1,2:Epm) = (uhalfx(1,1,2:Epm)-uhalfx(1,Enm,2:Epm)...
    +uhalfz(1,1,2:Epm)-uhalfz(1,1,1:Epmm))/dt/eust; %creases in wall
rhsfixed(1,En,2:Epm)=rhsfixed(1,1,2:Epm);

rhsfixed(Em,1,2:Epm)=(uhalfx(Em,1,2:Epm)-uhalfx(Em,Enm,2:Epm)...
    +uhalfz(Em,En,2:Epm)-uhalfz(Em,En,1:Epmm))/dt/eust; %top boundary
rhsfixed(Em,En,2:Epm)=rhsfixed(Em,1,2:Epm);

rhsfixed(1,1,1)=(uhalfx(1,1,1)-uhalfx(1,Enm,1)...
    +uhalfz(1,1,1)-uhalfz(1,1,Epm))/dt/eust;%creases at periodic boundary
rhsfixed(1,En,1)=rhsfixed(1,1,1);
rhsfixed(1,En,Ep)=rhsfixed(1,1,1);
rhsfixed(1,1,Ep)=rhsfixed(1,1,1);

rhsfixed(Em,1,1)=(uhalfx(Em,1,1)-uhalfx(Em,Enm,1)...
    +uhalfz(Em,En,1)-uhalfz(Em,En,Epm))/dt/eust; %top crease
rhsfixed(Em,En,1)=rhsfixed(Em,1,1);
rhsfixed(Em,En,Ep)=rhsfixed(Em,1,1);
rhsfixed(Em,1,Ep)=rhsfixed(Em,1,1);

rhsfixed(:,:,1)=rhsfixed(:,:,Ep);
rhsfixed(:,En,:)=rhsfixed(:,1,:);
rhsfixed=rhsfixed/h;

rhsfixed2=rhsfixed;
rhsfixed2(1,:,:)=1/2*rhsfixed2(1,:,:);
rhsfixed2(Em,:,:)=1/2*rhsfixed2(Em,:,:);
rhsfixed2=rhsfixed2-mean(mean(mean(rhsfixed2(:,1:Enm,1:Epm))));
    
psave=fmvPRESSUREprod3Dper4(p,rhsfixed2,h,Edens,levels,pcoef);

psave=psave-mean(mean(mean(psave(:,1:Enm,1:Epm))));

perr=max(max(max(abs(residualPRESSUREprod3Dper4(psave,rhsfixed2,pcoef{1})))));

% 
pnew=psave;
perrsave=perr;
count1=0;
perrkeep=perrsave;
c=0;
% pp(1)=perr;
while perr>allowerr && count1<1 && c<6
%     pp(c+1)=perr;
tic
    for c5=1:5
        pnew=mvPRESSUREprod3Dper4(pnew,rhsfixed2,h,Edens,levels,pcoef);
        pnew=pnew-mean(mean(mean(pnew(:,1:Enm,1:Epm))));
    end
    toc
    
    r=residualPRESSUREprod3Dper4(pnew,rhsfixed2,pcoef{1});
    perr=max(max(max(abs(r))));
    if perr<perrsave
        psave=pnew;
        perrsave=perr;
    else
    count1=count1+1;
    end
    perrkeep=[perrkeep;perrsave];
    c=c+1;
    
end
