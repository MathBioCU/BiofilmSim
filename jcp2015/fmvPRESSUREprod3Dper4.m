function phnew=fmvPRESSUREprod3Dper4(ph,fh,h,Edens,levels,pcoef)

lev=log2(h/(Edens{1}.dx))+1;
h2=2*h;
% [x,y]=meshgrid(0:dx:300,0:dx:300);
% Edens=transferLtoEdens(dx,x,y,X,initdensity,addldens,d0mean);

if lev==levels
    phnew=mvPRESSUREprod3Dper4(ph,fh,h,Edens,levels,pcoef);
else
%     [x2h,y2h]=meshgrid(0:dx2:300*10^-6,0:dx2:300*10^-6);
%     Edens2h=transferLtoEdens(dx2,x2h,y2h,X,initdensity,addldens,d0mean);
    f2h=restricthto2h3Dper4(residualPRESSUREprod3Dper2(ph,fh,pcoef{lev}));
    f2h=f2h-mean(mean(mean(f2h(:,1:end-1,1:end-1))));
    
    p2h=zeros((size(ph)+1)/2);
    p2h=fmvPRESSUREprod3Dper4(p2h,f2h,h2,Edens,levels,pcoef);
    phnew=ph+interpolate2htoh3D(p2h);
    phnew=phnew-mean(mean(mean(phnew(:,1:end-1,1:end-1))));
    phnew=mvPRESSUREprod3Dper4(phnew,fh,h,Edens,levels,pcoef);
end

