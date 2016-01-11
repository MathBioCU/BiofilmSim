function [fStrain,S,S0] = Calculate_Strains3(S,ux,uy,uz,S0,dt,st,xlength,ylength,zlength,h,x,y,z)

S1=S;
S1(:,3)=mod(S1(:,3),zlength);
S1(:,1)=mod(S1(:,1),xlength);
S1(S1(:,2)>ylength,2)=ylength;

Us=transferEtoLvel3Dper_e2a(h,ux,uy,uz,x,y,z,S1,zlength,xlength);
S=S+dt/st*Us;

V=S-S0; %displacement

% Sr=[S(:,1),S(:,2),S(:,3)+zlength];
% Sl=[S(:,1),S(:,2),S(:,3)-zlength];
% Vr=Sr-S0;
% Vl=Sl-S0;
% Sr=[S(:,1)+xlength,S(:,2),S(:,3)];
% Sl=[S(:,1)-xlength,S(:,2),S(:,3)];
% Vf=Sr-S0;
% Vb=Sl-S0;



LVz=zeros(100,1);
d1=LVz;
d2=LVz;
for i=1:100
    d1(i)=(V(i,3)-V(i+100,3))/(S(i,2)-S(i+100,2));
    d2(i)=(V(i+100,3)-V(i+200,3))/(S(i+100,2)-S(i+200,2));
end
for i=1:100
    if abs(d1(i)-d2(i))>0.1 && (abs(d1(i))>abs(mean(d1))+std(d1) || abs(d2(i))>abs(mean(d2))+std(d2))
        S1=S(i,:);
        if abs(d1)>abs(d2)
            S(i,:)=(S(i+100,:)-S(i+200,:))*(S(i,2)-S(i+100,2))/(S(i+100,2)-S(i+200,2))+S(i+100,:);
            S0(i,2)=S0(i+100,2)+(S0(i,2)-S0(i+100,2))*norm(S(i,:)-S(i+100,:))/norm(S1-S(i+100,:));
        elseif abs(d2)>abs(d1)
            S(i+200,:)=(S(i+100,:)-S(i,:))*(S(i+200,2)-S(i+100,2))/(S(i+100,2)-S(i,2))+S(i+100,:);
            S0(i+200,2)=S0(i+100,2)+(S0(i+200,2)-S0(i+100,2))*norm(S(i+200,:)-S(i+100,:))/norm(S1-S(i+100,:));
        end
    end
end
    
LVz=1/2*(d1+d2);
ind=find(LVz<mean(LVz)+std(LVz) & LVz>mean(LVz)-std(LVz));
LVz=LVz(ind);

fStrain=1/2*mean(LVz);


end

