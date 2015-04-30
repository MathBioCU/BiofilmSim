function [Efxnew,Efynew,Efznew]=transferLtoE3Dper_e2(h,x,y,z,X,sumFx,sumFy,sumFz,d0mean,zlength,xlength)

[Em,En,Ep]=size(x);
x=x(:,1:En-1,1:Ep-1);%this accounts for the fact that the first and last page are actually the same due to periodicity
y=y(:,1:En-1,1:Ep-1);
z=z(:,1:En-1,1:Ep-1); 


Efxnew=zeros(size(x));

Efynew=Efxnew;
Efznew=Efxnew;
rc=1/20;
coef=1/rc^3;
%change both to cubed for 3d

for c1=1:size(X,1)
    
        
    xi=ceil(X(c1,1)/h);
    yi=ceil(X(c1,2)/h);
    zi=mod(ceil(X(c1,3)/h),Ep);
    if zi==0 && X(c1,3)/h<Ep-1
        zi=Ep-1;
    end
    
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
    
%     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
    temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
    Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFx(c1)*coef*temp;
    Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFy(c1)*coef*temp;
    Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFz(c1)*coef*temp;
end

for c1=1:size(X,1)
    if X(c1,3)>(Ep-15)*h;
        
         xi=ceil(X(c1,1)/h);
        yi=ceil(X(c1,2)/h);
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
        zdist=abs(zlocal+zlength/h-X(c1,3))/rc;

    %     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
        temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
        
	Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFx(c1)*coef*temp;
        Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFy(c1)*coef*temp;
        Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFz(c1)*coef*temp;
    end
end

for c1=1:size(X,1)
    
        
    if X(c1,3)<h*15;
        
         xi=ceil(X(c1,1)/h);
        yi=ceil(X(c1,2)/h);
        zi=Ep;
        ximin=max(0,xi-10); ximax=min(En-2,xi+10); %change 15 to other values if dx changes (especially if it decreases)
        yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
        zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;


        xdist=abs(xlocal-X(c1,1))/rc;
        ydist=abs(ylocal-X(c1,2))/rc;
        zdist=abs(zlocal-zlength/h-X(c1,3))/rc;

    %     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
        temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
        Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFx(c1)*coef*temp;
        Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFy(c1)*coef*temp;
        Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFz(c1)*coef*temp;
    end
end

for c1=1:size(X,1)
    if X(c1,1)>(En-15)*h;
        
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

    %     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
        temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
        Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFx(c1)*coef*temp;
        Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFy(c1)*coef*temp;
        Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFz(c1)*coef*temp;
    end
end

for c1=1:size(X,1)
    
        
    if X(c1,1)<h*15;
        
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

    %     dtemp=sqrt(xdist.^2+ydist.^2+zdist.^2);
        temp=dirac_interp_mex(xdist).*dirac_interp_mex(ydist).*dirac_interp_mex(zdist);
        Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efxnew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFx(c1)*coef*temp;
        Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efynew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFy(c1)*coef*temp;
        Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=Efznew(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)+sumFz(c1)*coef*temp;
    end
end


Efxnew2=zeros(Em,En,Ep);
Efynew2=Efxnew2;
Efznew2=Efxnew2;
Efxnew2(:,1:En-1,1:Ep-1)=Efxnew; Efynew2(:,1:En-1,1:Ep-1)=Efynew; Efznew2(:,1:En-1,1:Ep-1)=Efznew;
Efxnew2(:,1:En-1,Ep)=Efxnew(:,1:En-1,1); Efynew2(:,1:En-1,Ep)=Efynew(:,1:En-1,1); Efznew2(:,1:En-1,Ep)=Efznew(:,1:En-1,1);
Efxnew2(:,En,1:Ep-1)=Efxnew(:,1,:); Efynew2(:,En,1:Ep-1)=Efynew(:,1,:); Efznew2(:,En,1:Ep-1)=Efznew(:,1,:);
Efxnew2(:,En,Ep)=Efxnew2(:,En,1); Efynew2(:,En,Ep)=Efynew2(:,En,1); Efznew2(:,En,Ep)=Efznew2(:,En,1);
Efxnew=Efxnew2; Efynew=Efynew2; Efznew=Efznew2;
