function Unew=transferEtoLvel3Dper_e2a(h,ux,uy,uz,x,y,z,X,zlength,xlength)
Unew=zeros(size(X));
[Em,En,Ep]=size(x);
x=x(:,1:En-1,1:Ep-1);%this accounts for the fact that the first and last page are actually the same due to periodicity
y=y(:,1:En-1,1:Ep-1);
z=z(:,1:En-1,1:Ep-1); 


rc=1/40;
coef=h^3/rc^3;

for c1=1:size(X,1)
    %form small cube around each position X(i) on the Eulerian grid.


    xi=mod(ceil(X(c1,1)/h),En-1);
    yi=mod(ceil(X(c1,2)/h),Em);
    zi=mod(ceil(X(c1,3)/h),Ep-1);
    
    if zi==0 && X(c1,3)/h<=Ep-1 && X(c1,3)>0
        zi=Ep-1;
    end
    if xi==0 && X(c1,1)/h<=En-1 && X(c1,1)>0
        xi=En-1;
    end
    
    if yi==0 && X(c1,2)/h>Em-1
        yi=Em-1;
    end
    
    ximin=max(0,xi-10); ximax=min(En-2,xi+10); %even though indices start at 1 and go to Ep-1, actual coordinates start at 0 and go to Ep-2 etc.
    yimin=max(0,yi-10); yimax=min(Em-1,yi+10);
    zimin=max(0,zi-10); zimax=min(Ep-2,zi+10);
    [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax, zimin:1:zimax);
    
    xdist=abs(h*xlocal-X(c1,1))/rc;
    ydist=abs(h*ylocal-X(c1,2))/rc;
    zdist=abs(h*zlocal-X(c1,3))/rc;  
        
    temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
    Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1).*temp))),...
        sum(sum(sum(uy(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1).*temp))),...
        sum(sum(sum(uz(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1).*temp)))];
    
        V=sum(sum(sum(temp)));
    V1=0; V2=0; V3=0; V4=0;
    if ximin>xi-10 
        [xlocal, ylocal, zlocal]=meshgrid(En-10:1:En-2,yimin:1:yimax,zimin:1:zimax);
        xdist=abs(h*xlocal-xlength-X(c1,1))/rc;
        ydist=abs(h*ylocal-X(c1,2))/rc;
        zdist=abs(h*zlocal-X(c1,3))/rc;  
        
        temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
        Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,En-9:En-1,zimin+1:zimax+1).*temp))),...
            sum(sum(sum(uy(yimin+1:yimax+1,En-9:En-1,zimin+1:zimax+1).*temp))),...
            sum(sum(sum(uz(yimin+1:yimax+1,En-9:En-1,zimin+1:zimax+1).*temp)))]+Unew(c1,:);
    
            V1=sum(sum(sum(temp)));
        
    end
     if ximax<xi+10 
        [xlocal, ylocal, zlocal]=meshgrid(0:1:9,yimin:1:yimax,zimin:1:zimax);
        xdist=abs(h*xlocal+xlength-X(c1,1))/rc;
        ydist=abs(h*ylocal-X(c1,2))/rc;
        zdist=abs(h*zlocal-X(c1,3))/rc;  
        
        temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
        Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,1:10,zimin+1:zimax+1).*temp))),...
            sum(sum(sum(uy(yimin+1:yimax+1,1:10,zimin+1:zimax+1).*temp))),...
            sum(sum(sum(uz(yimin+1:yimax+1,1:10,zimin+1:zimax+1).*temp)))]+Unew(c1,:);
    
            V2=sum(sum(sum(temp)));
        
     end
     if zimin>zi-10 
        [xlocal, ylocal, zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax,Ep-10:1:Ep-2);
        xdist=abs(h*xlocal-X(c1,1))/rc;
        ydist=abs(h*ylocal-X(c1,2))/rc;
        zdist=abs(h*zlocal-zlength-X(c1,3))/rc;  
        
        temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
        Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,ximin+1:ximax+1,Ep-9:Ep-1).*temp))),...
            sum(sum(sum(uy(yimin+1:yimax+1,ximin+1:ximax+1,Ep-9:Ep-1).*temp))),...
            sum(sum(sum(uz(yimin+1:yimax+1,ximin+1:ximax+1,Ep-9:Ep-1).*temp)))]+Unew(c1,:);
    
            V3=sum(sum(sum(temp)));
         if ximin>xi-10
            [xlocal, ylocal, zlocal]=meshgrid(En-10:1:En-2,yimin:1:yimax,Ep-10:1:Ep-2);
            xdist=abs(h*xlocal-xlength-X(c1,1))/rc;
            ydist=abs(h*ylocal-X(c1,2))/rc;
            zdist=abs(h*zlocal-zlength-X(c1,3))/rc;  

            temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
            Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,En-9:En-1,Ep-9:Ep-1).*temp))),...
                sum(sum(sum(uy(yimin+1:yimax+1,En-9:En-1,Ep-9:Ep-1).*temp))),...
                sum(sum(sum(uz(yimin+1:yimax+1,En-9:En-1,Ep-9:Ep-1).*temp)))]+Unew(c1,:);

            V3=sum(sum(sum(temp)))+V3;
        elseif ximax<xi+10  
            [xlocal, ylocal, zlocal]=meshgrid(0:1:9,yimin:1:yimax,Ep-10:1:Ep-2);
            xdist=abs(h*xlocal+xlength-X(c1,1))/rc;
            ydist=abs(h*ylocal-X(c1,2))/rc;
            zdist=abs(h*zlocal-zlength-X(c1,3))/rc;  

            temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
            Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,1:10,Ep-9:Ep-1).*temp))),...
                sum(sum(sum(uy(yimin+1:yimax+1,1:10,Ep-9:Ep-1).*temp))),...
                sum(sum(sum(uz(yimin+1:yimax+1,1:10,Ep-9:Ep-1).*temp)))]+Unew(c1,:);

            V3=sum(sum(sum(temp)))+V3;
        end
        
        
     end
     if zimax<zi+10
        [xlocal, ylocal, zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax,0:1:9);
        xdist=abs(h*xlocal-X(c1,1))/rc;
        ydist=abs(h*ylocal-X(c1,2))/rc;
        zdist=abs(h*zlocal+zlength-X(c1,3))/rc;  
        
        temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
        Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,ximin+1:ximax+1,1:10).*temp))),...
            sum(sum(sum(uy(yimin+1:yimax+1,ximin+1:ximax+1,1:10).*temp))),...
            sum(sum(sum(uz(yimin+1:yimax+1,ximin+1:ximax+1,1:10).*temp)))]+Unew(c1,:);
    
        V4=sum(sum(sum(temp)));
        if ximin>xi-10
            [xlocal, ylocal, zlocal]=meshgrid(En-10:1:En-2,yimin:1:yimax,0:1:9);
            xdist=abs(h*xlocal-xlength-X(c1,1))/rc;
            ydist=abs(h*ylocal-X(c1,2))/rc;
            zdist=abs(h*zlocal+zlength-X(c1,3))/rc;  

            temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
            Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,En-9:En-1,1:10).*temp))),...
                sum(sum(sum(uy(yimin+1:yimax+1,En-9:En-1,1:10).*temp))),...
                sum(sum(sum(uz(yimin+1:yimax+1,En-9:En-1,1:10).*temp)))]+Unew(c1,:);

            V4=sum(sum(sum(temp)))+V4;
        elseif ximax<xi+10  
            [xlocal, ylocal, zlocal]=meshgrid(0:1:9,yimin:1:yimax,0:1:9);
            xdist=abs(h*xlocal+xlength-X(c1,1))/rc;
            ydist=abs(h*ylocal-X(c1,2))/rc;
            zdist=abs(h*zlocal+zlength-X(c1,3))/rc;  

            temp=dirac_interp_new(xdist).*dirac_interp_new(ydist).*dirac_interp_new(zdist)*coef;
            Unew(c1,:)=[sum(sum(sum(ux(yimin+1:yimax+1,1:10,1:10).*temp))),...
                sum(sum(sum(uy(yimin+1:yimax+1,1:10,1:10).*temp))),...
                sum(sum(sum(uz(yimin+1:yimax+1,1:10,1:10).*temp)))]+Unew(c1,:);

            V4=sum(sum(sum(temp)))+V4;
        end
        
     end     
        Unew(c1,:)=Unew(c1,:)/(V+V1+V2+V3+V4);
       
end


