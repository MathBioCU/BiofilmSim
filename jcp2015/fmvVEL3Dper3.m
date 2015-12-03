function [vhx,vhy,vhz]=fmvVEL3Dper3(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef,v0)

lev=log2(h/(Edens{1}.dx))+1;
h2=2*h;
% [x,y]=meshgrid(0:dx:300,0:dx:300);
% Edens=transferLtoEdens(dx,x,y,X,initdensity,addldens,d0mean);
[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev},2); 
if lev==levels
    for i=1:1
    [vhx,vhy,vhz]=mvVEL3Dper2(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef);
    end
else
	fxres=zeros(size(vhx));
	fyres=zeros(size(vhx));
	fzres=zeros(size(vhx));
	[fxres,fyres,fzres]=residualvel3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev}, fxres,fyres,fzres);

    v2hx=zeros((size(vhx)+1)/2);
    v2hy=v2hx;
    v2hz=v2hx;
    v2hz(end,:,:)=v0;

    f2hx=restricthto2h3DVper2(fhx);
    f2hy=restricthto2h3DVper2(fhy);
    f2hz=restricthto2h3DVper2(fhz);

    [v2hx,v2hy,v2hz]=fmvVEL3Dper3(v2hx,v2hy,v2hz,f2hx,f2hy,f2hz,h2,Edens,levels,vcoef,v0);
     vhx=vhx+interpolate2htoh3D(v2hx);
     vhy=vhy+interpolate2htoh3D(v2hy);
     vhz=vhz+interpolate2htoh3D(v2hz);
     vhz(end,:,:)=v0;
      
    [vhx,vhy,vhz]=mvVEL3Dper2(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef);
%     for i=1:1
%     [vhx,vhy,vhz]=mvVEL3D(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef);
%     end
end

