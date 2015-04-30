clc
% clear all 
fr=[1 100];
%profile on

%for this test have to make sure transfer routines are updated so that
%correct region is chosen to add to. i.e. if need ximin=xi-10, ximax=xi+10;
%halving dx will required ximin=xi-20, ximax=xi+20 etc.
for i=1:2
    w=fr(i);
    dt=1/500;
    dx=1/64;
    dt=dt;
    numtimesteps=100;
    	
	fmax=100000;
    addlvisc=00;
    b=0;
    E=0;
    B01=b;
    E01=E; 
    sigma0=0.1;
    connectdist=3/18.5;
    levelsV=5;
    levelsP=5;
    
	JMain3DSim2v          

    vShear3=vShear;
    eShear3=eShear;
    fStrain2=fStrain;
    ux2=ux;
    uy2=uy;
    uz2=uz;
    p2=p;
    X2=X;
    A2=A; 
	Xs2=Xstore;
        
    b=B01;
	E=E01;
	dt=dt*2;
    numtimesteps=numtimesteps/2;
    
   
	JMain3DSim2v 

    ux5=ux;
    uy5=uy;
    uz5=uz;
    p5=p;
    X5=X;
	A5=A;      
	Xs5=Xstore;   
    vShear5=vShear;
    eShear5=eShear;
    fStrain5=fStrain;

    b=B01;
	E=E01;
	dt=dt/4;
    numtimesteps=numtimesteps*4;
    
    
	JMain3DSim2v 

    ux6=ux;
    uy6=uy;
    uz6=uz;
    p6=p;
    X6=X;
	A6=A;
	Xs6=Xstore;
    vShear6=vShear;
    eShear6=eShear;
    fStrain6=fStrain;

%         Eu2=sqrt(sum(sum(sum((restricthto2h3DVper2(uz6)-uz2).^2+(restricthto2h3DVper2(ux6)-ux2).^2+(restricthto2h3DVper2(uy6)-uy2).^2))))*h^(3/2);
%         Eu1=sqrt(sum(sum(sum((restricthto2h3DVper2(uz2)-uz5).^2+(restricthto2h3DVper2(ux2)-ux5).^2+(restricthto2h3DVper2(uy2)-uy5).^2))))*(2*h)^(3/2);
%         
        Ex2=sqrt(sum((X6(:,1)-X2(:,1)).^2+(X6(:,2)-X2(:,2)).^2+(X6(:,3)-X2(:,3)).^2));
        Ex1=sqrt(sum((X5(:,1)-X2(:,1)).^2+(X5(:,2)-X2(:,2)).^2+(X5(:,3)-X2(:,3)).^2));
        
        
       Eu2=sqrt(sum(sum(sum((uz6-uz2).^2+(ux6-ux2).^2+(uy6-uy2).^2))))*h^(3/2);
       Eu1=sqrt(sum(sum(sum((uz2-uz5).^2+(ux2-ux5).^2+(uy2-uy5).^2))))*h^(3/2);
       
       
       convvel=log2(Eu1/Eu2);
       convX=log2(Ex1/Ex2);

    	str1=num2str(w);
    	str2=num2str(fmax/1000);
    	str3=num2str(B01);
    	str4=num2str(1/(w*dt));
    	str5=num2str(connectdist);
    	str6=num2str(addlvisc);
    	str7=num2str(sigma0);
    	runid=['w',str1,'_f',str2,'_b',str3,'_cnd',str5,'_visc',str6,'_dt',str4,'_Verify_no_biofilm'];
    	save([runid,'.mat'],'uz2','ux2','uy2','X2','p2','ux5','uy5','uz5','p5','X5','ux6','uy6','uz6','p6','X6','h','Xstore','vShear','eShear', 'Xs2','Xs5','Xs6','A2','A5','A6')
    	Complete='Complete';
        fprintf('velocity: %f position: %f %s\n',convvel, convX,Complete);	

end
%p=profile('info');
%save profile_results p;
%profile off
runid

