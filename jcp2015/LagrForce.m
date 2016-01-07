function [XFe,YFe,ZFe,XFet,YFet,ZFet,Break,rowind,colind,matind,Intx,Inty,Intz,Dtemp]=LagrForce(X,xlength,ylength,zlength,matind,colind,rowind,Break,v0,U...
    ,Xrv,Yrv,Zrv,dt,K,d0,b,Xdist,Ydist,Zdist,A,XFe,YFe,ZFe,stuckt,stuckb,t,Intx,Inty,Intz,Dtempm)
    
    Zdistr=zeros(size(Zdist));
    Zdistl=Zdistr;
    Xdistr=Zdistr; 
    Xdistl=Zdistl;
    %calculate forces at the new positions
    
    Xdist(matind)=X(colind,1)-X(rowind,1);
    Xdistr(matind)=X(colind,1)-X(rowind,1)+xlength;
    Xdistl(matind)=X(colind,1)-X(rowind,1)-xlength;
    Ydist(matind)=X(colind,2)-X(rowind,2);
    Zdist(matind)=X(colind,3)-X(rowind,3);
    Zdistr(matind)=X(colind,3)-X(rowind,3)+zlength;
    Zdistl(matind)=X(colind,3)-X(rowind,3)-zlength;
    
    zind=find(abs(Zdist)>abs(Zdistr));
    Zdist(zind)=Zdistr(zind);
    zind=find(abs(Zdist)>abs(Zdistl));
    Zdist(zind)=Zdistl(zind);
    
    xind=find(abs(Xdist)>abs(Xdistr));
    Xdist(xind)=Xdistr(xind);
    xind=find(abs(Xdist)>abs(Xdistl));
    Xdist(xind)=Xdistl(xind);
    
    DSq=Xdist.^2+Ydist.^2+Zdist.^2;
    D=sqrt(DSq);
    
    breakind=find(D>2*d0);
    if ~isempty(breakind)

        A(breakind)=0;
        [rowind,colind]=find(A);
        matind=find(A);
        Break=Break+1;
    end
    Dtemp=(D-d0)./D;
    Fx=K.*Xdist.*Dtemp;
    Fy=K.*Ydist.*Dtemp;
    Fz=K.*Zdist.*Dtemp;
    
    Fx(~A)=0;
    Fy(~A)=0;
    Fz(~A)=0;
    Fx(isnan(Fx))=0;
    Fy(isnan(Fy))=0;
    Fz(isnan(Fz))=0;  
  
    Xrv(matind)=v0*(U(colind,1)-U(rowind,1));
    Yrv(matind)=v0*(U(colind,2)-U(rowind,2));
    Zrv(matind)=v0*(U(colind,3)-U(rowind,3));    
    aa=1;
    XrelVel=2*(1-(1./(1+abs(Xrv).^(1-aa))).*exp(-abs(Xrv).^(1-aa))).*abs(Xrv).^aa.*sign(Xrv);
    YrelVel=2*(1-(1./(1+abs(Yrv).^(1-aa))).*exp(-abs(Yrv).^(1-aa))).*abs(Yrv).^aa.*sign(Yrv);
    ZrelVel=2*(1-(1./(1+abs(Zrv).^(1-aa))).*exp(-abs(Zrv).^(1-aa))).*abs(Zrv).^aa.*sign(Zrv);

    VrelDotXrel=(XrelVel.*Xdist+YrelVel.*Ydist+ZrelVel.*Zdist)./DSq;
    Dx=(VrelDotXrel.*Xdist./D); Dx(isnan(Dx))=0; XFb=b.*Dx;
    Dy=(VrelDotXrel.*Ydist./D); Dy(isnan(Dy))=0; YFb=b.*Dy;
    Dz=(VrelDotXrel.*Zdist./D); Dz(isnan(Dz))=0; ZFb=b.*Dz;

    XFb(~A)=0;   XFb(isnan(XFb))=0;
    YFb(~A)=0;   YFb(isnan(YFb))=0;
    ZFb(~A)=0;   ZFb(isnan(ZFb))=0;
    
    dXFe=sparse(zeros(size(D)));
    dYFe=sparse(zeros(size(D)));
    dZFe=sparse(zeros(size(D)));
    dXFe(matind)=(XFe(colind)-XFe(rowind)).*(Xdist(matind))./D(matind)./d0(matind);
    dYFe(matind)=(YFe(colind)-YFe(rowind)).*(Ydist(matind))./D(matind)./d0(matind);
    dZFe(matind)=(ZFe(colind)-ZFe(rowind)).*(Zdist(matind))./D(matind)./d0(matind);

    Int=3;
    M=3;
    
    if M==1
        %Hooke's Law
        XFet=Fx; YFet=Fy; ZFet=Fz;
        Fx(stuckt,:)=0; Fy(stuckt,:)=0; Fz(stuckt,:)=0;
        Fx(stuckb,:)=0; Fy(stuckb,:)=0; Fz(stuckb,:)=0;
        XFe=Fx; YFe=Fy; ZFe=Fz;
        Intx=0; Inty=0; Intz=0;
    elseif M==2
        %Kelvin-Voigt Model
        XFet=Fx+XFb; YFet=Fy+YFb; ZFet=Fz+ZFb;
        Fx(stuckt,:)=0; Fy(stuckt,:)=0; Fz(stuckt,:)=0;
        Fx(stuckb,:)=0; Fy(stuckb,:)=0; Fz(stuckb,:)=0;
        XFb(stuckt,:)=0; YFb(stuckt,:)=0; ZFb(stuckt,:)=0;
        XFb(stuckb,:)=0; YFb(stuckb,:)=0; ZFb(stuckb,:)=0;
        XFe=Fx+XFb; YFe=Fy+YFb; ZFe=Fz+ZFb;
        Intx=0; Inty=0; Intz=0;
    elseif M==3
    
    %nonintegral form of Maxwell model
        if Int==1
            XFet=1./(b+K*dt).*(dt*b.*K.*Dx+b.*dXFe);
            YFet=1./(b+K*dt).*(dt*b.*K.*Dy+b.*dYFe);
            ZFet=1./(b+K*dt).*(dt*b.*K.*Dz+b.*dZFe);

            XFet(~A)=0;    Fx(stuckt,:)=0;     Fx(stuckb,:)=0;
            YFet(~A)=0;    Fy(stuckt,:)=0;     Fy(stuckb,:)=0;
            ZFet(~A)=0;    Fz(stuckt,:)=0;     Fz(stuckb,:)=0;
            XFet(isnan(XFet))=0;    XFb(stuckb,:)=0;    XFb(stuckt,:)=0;   
            YFet(isnan(YFet))=0;    YFb(stuckb,:)=0;    YFb(stuckt,:)=0;
            ZFet(isnan(ZFet))=0;    ZFb(stuckb,:)=0;    ZFb(stuckt,:)=0;
            dXFe(isnan(dXFe))=0;    dXFe(stuckb,:)=0;    dXFe(stuckt,:)=0;   
            dYFe(isnan(dYFe))=0;    dYFe(stuckb,:)=0;    dYFe(stuckt,:)=0;
            dZFe(isnan(dZFe))=0;    dZFe(stuckb,:)=0;    dZFe(stuckt,:)=0;
            Dx(isnan(Dx))=0;    Dx(stuckb,:)=0;    Dx(stuckt,:)=0;   
            Dy(isnan(Dy))=0;    Dy(stuckb,:)=0;    Dy(stuckt,:)=0;
            Dz(isnan(Dz))=0;    Dz(stuckb,:)=0;    Dz(stuckt,:)=0;

            XFe=1./(b+K*dt).*(dt*b.*K.*Dx+b.*dXFe);
            YFe=1./(b+K*dt).*(dt*b.*K.*Dy+b.*dYFe);
            ZFe=1./(b+K*dt).*(dt*b.*K.*Dz+b.*dZFe);

            XFe(~A)=0;  XFe(isnan(XFe))=0;
            YFe(~A)=0;  YFe(isnan(YFe))=0;  
            ZFe(~A)=0;  ZFe(isnan(ZFe))=0;
            Intx=0; Inty=0; Intz=0;
        elseif Int==2

            %integral form of maxwell model
            Intx=Intx+XFb.*exp(t*b./K)*dt;
            Inty=Inty+YFb.*exp(t*b./K)*dt;
            Intz=Intz+ZFb.*exp(t*b./K)*dt;

            XFet=exp(-t*b./K).*Intx;
            YFet=exp(-t*b./K).*Inty;
            ZFet=exp(-t*b./K).*Intz;

            XFet(~A)=0;    Fx(stuckt,:)=0;     Fx(stuckb,:)=0;
            YFet(~A)=0;    Fy(stuckt,:)=0;     Fy(stuckb,:)=0;
            ZFet(~A)=0;    Fz(stuckt,:)=0;     Fz(stuckb,:)=0;
            XFet(isnan(XFet))=0;    XFb(stuckb,:)=0;    XFb(stuckt,:)=0;   
            YFet(isnan(YFet))=0;    YFb(stuckb,:)=0;    YFb(stuckt,:)=0;
            ZFet(isnan(ZFet))=0;    ZFb(stuckb,:)=0;    ZFb(stuckt,:)=0;
            Intx(~A)=0;    Inty(~A)=0;         Intz(~A)=0;
            Intx(isnan(Intx))=0; Inty(isnan(Inty))=0; Intz(isnan(Intz))=0;

            XFe=exp(-t*b./K).*Intx;
            YFe=exp(-t*b./K).*Inty;
            ZFe=exp(-t*b./K).*Intz;

            XFe(~A)=0;  XFe(isnan(XFe))=0;
            YFe(~A)=0;  YFe(isnan(YFe))=0;  
            ZFe(~A)=0;  ZFe(isnan(ZFe))=0;
            XFe(stuckt,:)=0; XFe(stuckb,:)=0;
            YFe(stuckt,:)=0; YFe(stuckb,:)=0;
            ZFe(stuckt,:)=0; ZFe(stuckb,:)=0;
        elseif Int==3
            
            %integral form of maxwell model with integration by parts and
            %piecewise linear approximation of the strain 
            L=K./b;
            L(~A)=0;
            
            Int=Intx;
            Int=Int.*exp(-L*dt)+exp(-L*dt).*(Dtemp-Dtempm.*(1+dt*L)+exp(dt*L).*(Dtempm+Dtemp.*(-1+dt*L)))./(L.*L*dt);

            Int(~A)=0; 
            Int(isnan(Int))=0;
            
            c=2.0; %1.00 for Maxwell, >1 for Zener (1.05)
            XFet=(c*K.*Dtemp-K.*L.*Int).*Xdist;
            YFet=(c*K.*Dtemp-K.*L.*Int).*Ydist;
            ZFet=(c*K.*Dtemp-K.*L.*Int).*Zdist;
            
            XFet(~A)=0;  XFet(isnan(XFet))=0;
            YFet(~A)=0;  YFet(isnan(YFet))=0;  
            ZFet(~A)=0;  ZFet(isnan(ZFet))=0;
            
            XFe=XFet;
            YFe=YFet;
            ZFe=ZFet;
            XFe(stuckt,:)=0; XFe(stuckb,:)=0;
            YFe(stuckt,:)=0; YFe(stuckb,:)=0;
            ZFe(stuckt,:)=0; ZFe(stuckb,:)=0;
            Intx=Int;
        end
    end