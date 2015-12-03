function phnew=fmvPRESSUREprod3Dper2(ph,fh,h,Edens,levels,pcoef)

lev=log2(h/(Edens{1}.dx))+1;
h2=2*h;
% [x,y]=meshgrid(0:dx:300,0:dx:300);
% Edens=transferLtoEdens(dx,x,y,X,initdensity,addldens,d0mean);

fh(1,:,:)=1/2*fh(1,:,:);
fh(end,:,:)=1/2*fh(end,:,:);
fh=fh-mean(mean(mean(fh(:,1:end-1,1:end-1))));

if lev==levels
    phnew=mvPRESSUREprod3Dper2(ph,fh,h,Edens,levels,pcoef);
else
%     [x2h,y2h]=meshgrid(0:dx2:300*10^-6,0:dx2:300*10^-6);
%     Edens2h=transferLtoEdens(dx2,x2h,y2h,X,initdensity,addldens,d0mean);
    f2h=restricthto2h3Dper2(residualPRESSUREprod3Dper2(ph,fh,pcoef{lev}));
    
    p2h=zeros((size(ph)+1)/2);
    p2h=fmvPRESSUREprod3Dper2(p2h,f2h,h2,Edens,levels,pcoef);
    phnew=ph+interpolate2htoh3D(p2h);
    phnew=phnew-mean(mean(mean(phnew(:,1:end-1,1:end-1))));
    phnew=mvPRESSUREprod3Dper2(phnew,fh,h,Edens,levels,pcoef);
end

