function [vhx,vhy,vhz]=fmvVEL3Dper4(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef)

lev=log2(h/(Edens{1}.dx))+1;
h2=2*h;
% [x,y]=meshgrid(0:dx:300,0:dx:300);
% Edens=transferLtoEdens(dx,x,y,X,initdensity,addldens,d0mean);
if lev==1
	[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{1},2);
else
	[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev},2); 
end
if lev==levels
    for i=1:1
    [vhx,vhy,vhz]=mvVEL3Dper4(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef);
    end
else
if lev==1
	fxres=zeros(size(vhx));
	fyres=zeros(size(vhx));
	fzres=zeros(size(vhx));
	[fxres,fyres,fzres]=residualvel3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{1}, fxres,fyres,fzres);
else
	fxres=zeros(size(vhx));
	fyres=zeros(size(vhx));
	fzres=zeros(size(vhx));
	[fxres,fyres,fzres]=residualvel3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev}, fxres,fyres,fzres);
end

    v2hx=zeros((size(vhx)+1)/2);
    v2hy=v2hx;
    v2hz=v2hx;
%     f2hxres=restricthto2h3DV(fxres);
%     f2hyres=restricthto2h3DV(fyres);
%     f2hzres=restricthto2h3DV(fzres);
%     for i=1:1
%     [v2hx,v2hy,v2hz]=fmvVEL3D(v2hx,v2hy,v2hz,f2hxres,f2hyres,f2hzres,h2,Edens,levels,vcoef);
%     end
%     [vhx,vhy,vhz]=relaxVEL3D(vhx+interpolate2htoh3D(v2hx),vhy+interpolate2htoh3D(v2hy),vhz+interpolate2htoh3D(v2hz),fhx,fhy,fhz,2,vcoef{lev});
    [v2hx,v2hy,v2hz]=fmvVEL3Dper4(v2hx,v2hy,v2hz,restricthto2h3DVper2_mex(fxres),restricthto2h3DVper2_mex(fyres),restricthto2h3DVper2_mex(fzres),h2,Edens,levels,vcoef);
    [vhx,vhy,vhz]=mvVEL3Dper4(vhx+interpolate2htoh3D(v2hx),vhy+interpolate2htoh3D(v2hy),vhz+interpolate2htoh3D(v2hz),fhx,fhy,fhz,h,Edens,levels,vcoef);
%     for i=1:1
%     [vhx,vhy,vhz]=mvVEL3D(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef);
%     end
end

