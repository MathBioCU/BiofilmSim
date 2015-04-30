function phnew=mvPRESSUREprod3Dper2(ph,fh,h,Edens,levels,pcoef)

lev=log2(h/(Edens{1}.dx))+1;
h2=2*h;
numrel1=2;
numrel2=2;
% [x,y]=meshgrid(0:dx:300,0:dx:300);
% Edens=transferLtoEdens(dx,x,y,X,initdensity,addldens,d0mean);
ph=relaxPRESSUREprodRB3Dper2(ph,fh,numrel1,pcoef{lev});
if lev==levels%30
    phnew=relaxPRESSUREprodRB3Dper2(ph,fh,10,pcoef{lev});
else
    f2h=restricthto2h3Dper2(residualPRESSUREprod3Dper2(ph,fh,pcoef{lev}));
    f2h(1,:,:)=1/2*f2h(1,:,:);
    f2h(end,:,:)=1/2*f2h(end,:,:);
    f2h=f2h-mean(mean(mean(f2h(:,1:end-1,1:end-1))));

%     oneh=ones(size(f2h));
%     oneh(1,:,:)=1/2;
%     oneh(end,:,:)=1/2;
%     f2h=f2h-sum(sum(sum(f2h.*oneh)))/sum(sum(sum(oneh.*oneh)))*oneh;

    p2h=zeros((size(ph)+1)/2);
%     [x,y]=meshgrid(0:dx2:300*10^-6,0:dx2:300*10^-6);
%     Edens2h=transferLtoEdens(dx2,x,y,X,initdensity,addldens,d0mean);
    p2h=mvPRESSUREprod3Dper2(p2h,f2h,h2,Edens,levels,pcoef);
    phnew=ph+interpolate2htoh3D(p2h);
    phnew=relaxPRESSUREprodRB3Dper2(phnew,fh,numrel2,pcoef{lev});
    phnew=phnew-mean(mean(mean(phnew(:,1:end-1,1:end-1))));
end
