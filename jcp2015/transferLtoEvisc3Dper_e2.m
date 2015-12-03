function viscmat=transferLtoEvisc3Dper_e2(h,x,y,z,X,visc,addlvisc,numOfnonzero,lev,zlength,xlength)


[Em,En,Ep]=size(x);
x=x(:,1:En-1,1:Ep-1);%this accounts for the fact that the first and last page are actually the same due to periodicity
y=y(:,1:En-1,1:Ep-1);
z=z(:,1:En-1,1:Ep-1); 

viscmat=visc*ones(size(x));


rc=1/20;
rc=rc*2^(lev-1);
% dzero=.1;%this is the distance that the viscosity covers (the radius of the support of viscosity) .1=5microns 

coefvisc=addlvisc*8;%times 4 b/c dirac_interp_mex(0)=.5 so .5^3=1/8 and we want max to be 500
if lev>1
    coefvisc=coefvisc/(2^(lev-1))^2;
end
% coef1=addlvisc;
for c1=1:numOfnonzero
          
   xi=mod(ceil(X(c1,1)/h),En);
    yi=mod(ceil(X(c1,2)/h),Em);
    zi=mod(ceil(X(c1,3)/h),Ep);
    
    ximin=max(0,xi-10); ximax=min(En-2,xi+10); %change 15 to other values if dx changes (especially if it decreases)
    yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
    zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal-X(c1,3))/rc;
       
    
if addlvisc>0
    if lev<=10% maxvisc>(visc+addlvisc)
        
tempval=1*coefvisc*(dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist));
viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=tempval+viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);

    end




end
end

for c1=1:numOfnonzero
    if X(c1,3)>(Ep-20)*h
          
    xi=mod(ceil(X(c1,1)/h),En);
    yi=mod(ceil(X(c1,2)/h),Em);
    zi=0;
    
    ximin=max(0,xi-10); ximax=min(En-2,xi+10); %change 15 to other values if dx changes (especially if it decreases)
    yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
    zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal+zlength-X(c1,3))/rc;
    
    if addlvisc>0
        if lev<=10% maxvisc>(visc+addlvisc)
        tempval=1*coefvisc*(dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist));

            viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=tempval+viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);

        end




    end
    end
end

for c1=1:numOfnonzero
    if X(c1,3)<20*h
             
        xi=mod(ceil(X(c1,1)/h),En);
        yi=mod(ceil(X(c1,2)/h),Em);
        zi=Ep-1;

        ximin=max(0,xi-10); ximax=min(En-2,xi+10); %change 15 to other values if dx changes (especially if it decreases)
        yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
        zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;


        xdist=abs(xlocal-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-zlength-X(c1,3))/rc;

        if addlvisc>0
            if lev<=10% maxvisc>(visc+addlvisc)
            tempval=1*coefvisc*(dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist));

                viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=tempval+viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);

            end
        end
    end

end

for c1=1:numOfnonzero
    if X(c1,1)>(En-20)*h
          
    xi=0;
    yi=mod(ceil(X(c1,2)/h),Em);
    zi=mod(ceil(X(c1,3)/h),Ep);
    
    ximin=max(0,xi-10); ximax=min(En-2,xi+10); %change 15 to other values if dx changes (especially if it decreases)
    yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
    zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal+xlength-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal-X(c1,3))/rc;
    
    if addlvisc>0
        if lev<=10% maxvisc>(visc+addlvisc)
        tempval=1*coefvisc*(dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist));

            viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=tempval+viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);

        end




    end
    end
end

for c1=1:numOfnonzero
    if X(c1,1)<20*h
             
        xi=En;
        yi=mod(ceil(X(c1,2)/h),Em);
        zi=mod(ceil(X(c1,3)/h),Ep);

        ximin=max(0,xi-10); ximax=min(En-2,xi+10); %change 15 to other values if dx changes (especially if it decreases)
        yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
        zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;


        xdist=abs(xlocal-xlength-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-X(c1,3))/rc;

        if addlvisc>0
            if lev<=10% maxvisc>(visc+addlvisc)
            tempval=1*coefvisc*(dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist));

                viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=tempval+viscmat(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);

            end
        end
    end

end






viscmat2=zeros(Em,En,Ep);
viscmat2(:,1:En-1,1:Ep-1)=viscmat;
viscmat2(:,1:En-1,Ep)=viscmat(:,1:En-1,1);
viscmat2(:,En,:)=viscmat2(:,1,:);
viscmat=viscmat2;


