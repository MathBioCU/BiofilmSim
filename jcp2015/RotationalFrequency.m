%compute average rotational frequency of bacteria in aggregate is at
%tumbles

mX=zeros(c3+1,3);
for i=1:c3+1
    for j=0:2
        mX(i,j+1)=mean(Xrstore(:,3*i-2+j));
    end
end

Xc=zeros(c3+1,3);
rot_freq=zeros(length(X),1);
for i=1:length(X)
    for j=1:c3+1
        for k=0:2
            Xc(j,k+1)=Xrstore(i,3*j-2+k)-mX(j,k+1);
        end
    end
    Xc(:,3)=Xc(:,3)-mean(Xc(:,3));
    t0=zeros(20,1);
    
    k=2;
    kk=1;
    while k<c3
        if Xc(k,3)>0 && Xc(k-1,3)<0 
            t0(kk)=k*dt;
            kk=kk+1;
        end
        k=k+1; 
    end
    
    t1=zeros(19,1);
    for j=1:19
        t1(j)=t0(j+1)-t0(j);
    end
    t1=t1(t1>0);
    rot_freq(i)=mean(t1);
    
end
    
    
    