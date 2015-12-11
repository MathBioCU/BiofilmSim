function [vhx,vhy,vhz]=mvVEL3Dper4(vhx,vhy,vhz,fhx,fhy,fhz,h,Edens,levels,vcoef)

lev=log2(h/(Edens{1}.dx))+1;

h2=2*h;
numrel1=1;
numrel2=2;
if lev==1
	[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{1},numrel1);	
else
	[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev},numrel1);
end
if lev==levels
    [vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev},20);
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
    	[fxres,fyres,fzres]=residualvel3Dper3_mex(vhx,vhy,vhz,fhx,fhy,fhz,vcoef{lev},fxres,fyres,fzres);
    end
    v2hx=zeros((size(vhx)+1)/2);
    v2hy=v2hx;
    v2hz=v2hx;

    [v2hx,v2hy,v2hz]=mvVEL3Dper4(v2hx,v2hy,v2hz,restricthto2h3DVper2(fxres),restricthto2h3DVper2(fyres),restricthto2h3DVper2(fzres),h2,Edens,levels,vcoef);
    if lev==1
    	[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx+interpolate2htoh3D(v2hx),vhy+interpolate2htoh3D(v2hy),vhz+interpolate2htoh3D(v2hz),fhx,fhy,fhz,vcoef{lev},numrel2);
    else
    	[vhx,vhy,vhz]=relaxVEL3Dper3_mex(vhx+interpolate2htoh3D(v2hx),vhy+interpolate2htoh3D(v2hy),vhz+interpolate2htoh3D(v2hz),fhx,fhy,fhz,vcoef{lev},numrel2);
    end
end

end
