function dru=div_rho_u(x,y,z,ux,uy,uz,X,U,addldens,h, xlength, ylength, zlength)
[Em,En,Ep]=size(x);
x=x(:,1:En-1,1:Ep-1);%this accounts for the fact that the first and last page are actually the same due to periodicity
y=y(:,1:En-1,1:Ep-1);
z=z(:,1:En-1,1:Ep-1); 

dru=zeros(size(x));

rc=1/30;
coef=1/rc^3;

for c1=1:length(X(:,1))
          
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


    xdist=(xlocal-X(c1,1))/rc;
    ydist=(ylocal-X(c1,2))/rc;
    zdist=(zlocal-X(c1,3))/rc;
       
    
    tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
    tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
    tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
    
    chi=ones(size(tempvalx));
    
    dru(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)).*tempvalx+...
        (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)).*tempvaly+...
        (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1)).*tempvalz+...
        dru(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:zimax+1);

    
    if ximin==0 
        ximin=round(En-15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,zimin:1:zimax);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;
        
        xdist=(xlocal-X(c1,1))/rc;
        ydist=(ylocal-X(c1,2))/rc;
        zdist=(zlocal-X(c1,3))/rc;


        tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
        tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
        tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
        
        chi=ones(size(tempvalx));

        dru(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:zimax+1)).*tempvalx+...
            (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:zimax+1)).*tempvaly+...
            (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:zimax+1)).*tempvalz+...
            dru(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:zimax+1);

               
        if zimin==0
            zimin=round(Ep-15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,zimin:1:Ep-2);
            
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
           xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
            
            chi=ones(size(tempvalx));

            dru(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)).*tempvalz+...
                dru(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1);
            
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
                        
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
           
            chi=ones(size(tempvalx));

            dru(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvalz+...
                dru(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1);

        elseif zimax==Ep-2
            zimax=round(15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,0:1:zimax);
           
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
             xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);

            chi=ones(size(tempvalx));


            dru(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)).*tempvalz+...
                dru(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1);
            
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,0:1:zimax);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
                        
            chi=ones(size(tempvalx));
            
            dru(yimin+1:yimax+1,1:ximax+1,1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,1:ximax+1,1:zimax+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,1:ximax+1,1:zimax+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,1:ximax+1,1:zimax+1)).*tempvalz+...
                dru(yimin+1:yimax+1,1:ximax+1,1:zimax+1);
        end
    elseif ximax==En-2
        ximax=round(15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;

         xdist=(xlocal-X(c1,1))/rc;
        ydist=(ylocal-X(c1,2))/rc;
        zdist=(zlocal-X(c1,3))/rc;


        tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
        tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
        tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);

            chi=ones(size(tempvalx));
            
            dru(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvalz+...
                dru(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1);
       
        if zimin==0
            zimin=round(Ep-15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
            
            chi=ones(size(tempvalx));

            dru(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1)).*tempvalz+...
                dru(yimin+1:yimax+1,1:ximax+1,zimin+1:Ep-2+1);
            
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,zimin:1:Ep-2);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
          
            chi=ones(size(tempvalx));

            dru(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1)).*tempvalz+...
                dru(yimin+1:yimax+1,ximin+1:En-2+1,zimin+1:Ep-2+1);
             
        elseif zimax==Ep-2
            zimax=round(15/(h/0.0137));
            [xlocal,ylocal,zlocal]=meshgrid(ximin:1:En-2,yimin:1:yimax,0:1:zimax);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
            
            chi=ones(size(tempvalx));

            dru(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1)).*tempvalz+...
                dru(yimin+1:yimax+1,ximin+1:En-2+1,1:zimax+1);
            
            [xlocal,ylocal,zlocal]=meshgrid(0:1:ximax,yimin:1:yimax,0:1:zimax);
            xlocal=h*xlocal;
            ylocal=h*ylocal;
            zlocal=h*zlocal;
        
            xdist=(xlocal-X(c1,1))/rc;
            ydist=(ylocal-X(c1,2))/rc;
            zdist=(zlocal-X(c1,3))/rc;


            tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
            tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
            tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);

            chi=ones(size(tempvalx));
            
            dru(yimin+1:yimax+1,1:ximax+1,1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,1:ximax+1,1:zimax+1)).*tempvalx+...
                (U(c1,2)*chi-uy(yimin+1:yimax+1,1:ximax+1,1:zimax+1)).*tempvaly+...
                (U(c1,3)*chi-uz(yimin+1:yimax+1,1:ximax+1,1:zimax+1)).*tempvalz+...
                dru(yimin+1:yimax+1,1:ximax+1,1:zimax+1);
        end
        
    elseif zimin==0
        zimin=round(Ep-15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax,zimin:1:Ep-2);
        xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;
        
          xdist=(xlocal-X(c1,1))/rc;
        ydist=(ylocal-X(c1,2))/rc;
        zdist=(zlocal-X(c1,3))/rc;


        tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
        tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
        tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);

        chi=ones(size(tempvalx));

        dru(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-2+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-2+1)).*tempvalx+...
            (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-2+1)).*tempvaly+...
            (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-2+1)).*tempvalz+...
            dru(yimin+1:yimax+1,ximin+1:ximax+1,zimin+1:Ep-2+1);
        
    elseif zimax==Ep-2
        zimax=round(15/(h/0.0137));
        [xlocal,ylocal,zlocal]=meshgrid(ximin:1:ximax,yimin:1:yimax,0:1:zimax);
       xlocal=h*xlocal;
        ylocal=h*ylocal;
        zlocal=h*zlocal;
    
        xdist=(xlocal-X(c1,1))/rc;
        ydist=(ylocal-X(c1,2))/rc;
        zdist=(zlocal-X(c1,3))/rc;


        tempvalx=dirac_interp_dr(xdist).*dirac_interp1(abs(ydist)).*dirac_interp1(abs(zdist));
        tempvaly=dirac_interp1(abs(xdist)).*dirac_interp_dr(ydist).*dirac_interp1(abs(zdist));
        tempvalz=dirac_interp1(abs(xdist)).*dirac_interp1(abs(ydist)).*dirac_interp_dr(zdist);
        
        chi=ones(size(tempvalx));

        dru(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1)=(U(c1,1)*chi-ux(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1)).*tempvalx+...
            (U(c1,2)*chi-uy(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1)).*tempvaly+...
            (U(c1,3)*chi-uz(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1)).*tempvalz+...
            dru(yimin+1:yimax+1,ximin+1:ximax+1,1:zimax+1);

    end  
    
end


drunew2=zeros(Em,En,Ep);
drunew2(:,1:En-1,1:Ep-1)=dru;
drunew2(:,1:En-1,Ep)=dru(:,1:En-1,1);
drunew2(:,En,1:Ep-1)=dru(:,1,:);
drunew2(:,En,Ep)=drunew2(:,En,1);
dru=addldens*8*drunew2;
