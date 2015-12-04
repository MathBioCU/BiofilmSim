function phnew=mvPRESSUREprod3Dper4(ph,fh,h,Edens,levels,pcoef)

lev=log2(h/(Edens{1}.dx))+1;
h2=2*h;
numrel1=2;
numrel2=2;

ph=ph-mean(mean(mean(ph(:,1:end-1,1:end-1))));
ph=relaxPRESSUREprodRB3Dper4(ph,fh,numrel1,pcoef{lev});
if lev==levels%30
     ph=ph-mean(mean(mean(ph(:,1:end-1,1:end-1))));
     phnew=relaxPRESSUREprodRB3Dper4(ph,fh,10,pcoef{lev});
else
    
    f2h=restricthto2h3Dper4(residualPRESSUREprod3Dper4(ph,fh,pcoef{lev}));
    f2h=f2h-mean(mean(mean(f2h(:,1:end-1,1:end-1))));

    p2h=zeros((size(ph)+1)/2);
    p2h=mvPRESSUREprod3Dper4(p2h,f2h,h2,Edens,levels,pcoef);
    phnew=ph+interpolate2htoh3D(p2h);
 
    phnew=phnew-mean(mean(mean(phnew(:,1:end-1,1:end-1))));
    phnew=relaxPRESSUREprodRB3Dper4(phnew,fh,numrel2,pcoef{lev});
 
end
