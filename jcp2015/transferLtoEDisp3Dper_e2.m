function Estrain=transferLtoEDisp3Dper_e2(h,x,y,z,X,U,numOfnonzero,zlength,xlength, rc)


[Em,En,Ep]=size(x);
x=x(:,1:En-1,1:Ep-1);%this accounts for the fact that the first and last page are actually the same due to periodicity
y=y(:,1:En-1,1:Ep-1);
z=z(:,1:En-1,1:Ep-1); 

Estrain=zeros(size(x));


%rc=1/20;
coef=1/rc^3;%change both to cubed for 3d


for c1=1:numOfnonzero
    
    xi=mod(ceil(X(c1,1)/h),En);
    yi=ceil(X(c1,2)/h);
    zi=mod(ceil(X(c1,3)/h),Ep);
    if zi==0 && X(c1,3)/h<Ep-1
        zi=Ep-1;
    end
    if xi==0 && X(c1,1)/h<En-1
        xi=En-1;
    end
    
    ximin=max(0,xi-15); ximax=min(En-2,xi+15); %change 15 to other values if dx changes (especially if it decreases)
    yimin=max(0,yi-15); yimax=min(Em-1,yi+15);
    zimin=max(0,zi-15); zimax=min(Ep-2,zi+15);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal-X(c1,3))/rc;
           
    Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+U(c1)*coef*dirac_interp_mex(xdist).*dirac_interp_mex(ydist)...
        .*dirac_interp_mex(zdist);
    
    
end
for c1=1:numOfnonzero         
    if X(c1,3)>(Ep-15)*h

    xi=mod(ceil(X(c1,1)/h),En);
    yi=ceil(X(c1,2)/h);
    zi=0;
    
    if xi==0 && X(c1,1)/h<En-1
        xi=En-1;
    end
    
    ximin=max(0,xi-20); ximax=min(En-2,xi+20);
    yimin=max(0,yi-20); yimax=min(Em-1,yi+20);
    zimin=max(0,zi-20); zimax=min(Ep-2,zi+20);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal+zlength-X(c1,3))/rc;
           
     Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+U(c1)*coef*dirac_interp_mex(xdist).*dirac_interp_mex(ydist)...
        .*dirac_interp_mex(zdist);
    
    end 
end
for c1=1:numOfnonzero
    if X(c1,3)<15*h
          
    xi=mod(ceil(X(c1,1)/h),En);
    yi=ceil(X(c1,2)/h);
    zi=Ep;
    if xi==0 && X(c1,1)/h<En-1
        xi=En-1;
    end
    
    ximin=max(0,xi-15); ximax=min(En-2,xi+15);
    yimin=max(0,yi-15); yimax=min(Em-1,yi+15);
    zimin=max(0,zi-15); zimax=min(Ep-2,zi+15);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    xlocal=h*xlocal;
    ylocal=h*ylocal;
    zlocal=h*zlocal;


    xdist=abs(xlocal-X(c1,1))/rc;
    ydist=abs(ylocal-X(c1,2))/rc;
    zdist=abs(zlocal-zlength-X(c1,3))/rc;
           
      Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+U(c1)*coef*dirac_interp_mex(xdist).*dirac_interp_mex(ydist)...
        .*dirac_interp_mex(zdist);
    
    
    end
end

for c1=1:size(X,1)
    if X(c1,1)>(En-15)*h;
        
        xi=0;
        yi=mod(ceil(X(c1,2)/h),Em);
        zi=mod(ceil(X(c1,3)/h),Ep);
        ximin=max(0,xi-15); ximax=min(En-2,xi+15); %change 15 to other values if dx changes (especially if it decreases)
        yimin=max(0,yi-15); yimax=min(Em-1,yi+15);
        zimin=max(0,zi-15); zimax=min(Ep-2,zi+15);
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;


        xdist=abs(xlocal+xlength-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-X(c1,3))/rc;

    %     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
        temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
        temp(xdist>2 | ydist>2 | zdist>2)=0;
       Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+U(c1)*coef*temp;
      
    end
end

for c1=1:size(X,1)
    
        
    if X(c1,1)<h*15;
        
        xi=En;
        yi=mod(ceil(X(c1,2)/h),Em);
        zi=mod(ceil(X(c1,3)/h),Ep);
        ximin=max(0,xi-15); ximax=min(En-2,xi+15); %change 15 to other values if dx changes (especially if it decreases)
        yimin=max(0,yi-15); yimax=min(Em-1,yi+15);
        zimin=max(0,zi-15); zimax=min(Ep-2,zi+15);
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;


        xdist=abs(xlocal-xlength-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-X(c1,3))/rc;

    %     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
        temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
        temp(xdist>2 | ydist>2 | zdist>2)=0;
        Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Estrain(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+U(c1)*coef*temp;
       
    end
end


Estrain2=zeros(Em,En,Ep);
Estrain2(:,1:En-1,1:Ep-1)=Estrain;
Estrain2(:,1:En-1,Ep)=Estrain(:,:,1);
Estrain2(:,1:En-1,1)=Estrain2(:,1:En-1,Ep);
Estrain2(:,En,1:Ep-1)=Estrain(:,1,:);
Estrain2(:,1,1:Ep-1)=Estrain2(:,En,1:Ep-1);
Estrain2(:,En,Ep)=Estrain2(:,1,Ep);
Estrain=Estrain2;
