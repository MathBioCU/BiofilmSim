function Edens=transferLtoEdens3Dper_e2a(h,x,y,z,X,initdensity,addldens,d0mean,numOfnonzero,zlength,xlength)



[Em,En,Ep]=size(x);
x=x(:,1:En-1,1:Ep-1);%this accounts for the fact that the first and last page are actually the same due to periodicity
y=y(:,1:En-1,1:Ep-1);
z=z(:,1:En-1,1:Ep-1); 

Edens=initdensity*ones(size(x));


rc=1/30;
coef=addldens*8;%change both to cubed for 3d;

% coef1=addldens;
for c1=1:numOfnonzero
          
    xi=mod(ceil(X(c1,1)/h),En);
    yi=mod(ceil(X(c1,2)/h),Em);
    zi=mod(ceil(X(c1,3)/h),Ep);
    
    ximin=max(0,floor(xi-10/(h/0.0137))); ximax=min(En-2,ceil(xi+10/(h/0.0137))); %change 15 to other values if dx changes (especially if it decreases)
    yimin=max(0,floor(yi-10/(h/0.0137))); yimax=min(Em-1,ceil(yi+10/(h/0.0137)));
    zimin=max(0,floor(zi-10/(h/0.0137))); zimax=min(Ep-2,ceil(zi+10/(h/0.0137)));
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal-X(c1,3))/rc;
           
    tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
    Edens(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=tempval+Edens(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);
    
    if ximin==0 
        ximin=round(En-15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;
        
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-X(c1,3))/rc;
        xdist=abs(xlocal-X(c1,1)-xlength)/rc;
        
        tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
        Edens(yimin+1:yimax+1,ximin+1:En-1,zimin+1:zimax+1)=tempval+Edens(yimin+1:yimax+1,ximin+1:En-1,zimin+1:zimax+1);
       
        if zimin==0
            zimin=round(Ep-15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,zimin:1:Ep-2);
            
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1)-xlength)/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)-zlength)/rc;
            
            tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
            Edens(yimin+1:yimax+1,ximin+1:En-1,zimin+1:Ep-1)=tempval+Edens(yimin+1:yimax+1,ximin+1:En-1,zimin+1:Ep-1);
            
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
                        
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1))/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)-zlength)/rc;
            
             tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
             Edens(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-1)=tempval+Edens(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-1);
        elseif zimax==Ep-2
            zimax=round(15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,0:1:zimax);
           
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1)-xlength)/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)+zlength)/rc;
            
            tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
            Edens(yimin+1:yimax+1,ximin+1:En-1,1:zimax+1)=tempval+Edens(yimin+1:yimax+1,ximin+1:En-1,1:zimax+1);
            
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,0:1:zimax);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1))/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)+zlength)/rc;
            
             tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
             Edens(yimin+1:yimax+1,1:ximax+1,1:zimax+1)=tempval+Edens(yimin+1:yimax+1,1:ximax+1,1:zimax+1);
        end
    elseif ximax==En-2
        ximax=round(15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,zimin:1:zimax);
        xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1)+xlength)/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3))/rc;
        
        tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
        Edens(yimin+1:yimax+1,1:1:ximax+1,zimin+1:zimax+1)=tempval+Edens(yimin+1:yimax+1,1:ximax+1,zimin+1:zimax+1);
       
        if zimin==0
            zimin=round(Ep-15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1)+xlength)/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)-zlength)/rc;
            
            tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
            Edens(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-1)=tempval+Edens(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-1);
            
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,zimin:1:Ep-2);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1))/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)-zlength)/rc;
            
             tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
             Edens(yimin+1:yimax+1,ximin+1:En-1,zimin+1:Ep-1)=tempval+Edens(yimin+1:yimax+1,ximin+1:En-1,zimin+1:Ep-1);
        elseif zimax==Ep-2
            zimax=round(15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,0:1:zimax);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1))/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)+zlength)/rc;
            
            tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
            Edens(yimin+1:yimax+1,ximin+1:En-1,1:zimax+1)=tempval+Edens(yimin+1:yimax+1,ximin+1:En-1,1:zimax+1);
            
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,0:1:zimax);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=abs(xlocal-X(c1,1)+xlength)/rc;
            ydist=abs(ylocal-X(c1,2))/rc;
            zdist=abs(zlocal-X(c1,3)+zlength)/rc;
            
             tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
             Edens(yimin+1:yimax+1,1:ximax+1,1:zimax+1)=tempval+Edens(yimin+1:yimax+1,1:ximax+1,1:zimax+1);
        end
        
    elseif zimin==0
        zimin=round(Ep-15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;
        
        xdist=abs(xlocal-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-X(c1,3)-zlength)/rc;
        
        tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
        Edens(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-1)=tempval+Edens(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-1);
        
    elseif zimax==Ep-2
        zimax=round(15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax,0:1:zimax);
       xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;
    
        xdist=abs(xlocal-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-X(c1,3)+zlength)/rc;
        
        tempval=coef*(dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist));
        Edens(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1)=tempval+Edens(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1);
    end  
    
end

Edens(Edens>addldens+initdensity)=addldens+initdensity;

Edens2=zeros(Em,En,Ep);
Edens2(:,1:En-1,1:Ep-1)=Edens;
Edens2(:,1:En-1,Ep)=Edens(:,1:En-1,1);
Edens2(:,En,:)=Edens2(:,1,:);
Edens=Edens2;


