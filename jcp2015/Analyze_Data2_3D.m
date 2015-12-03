function [G1, G2, Delta, maxshear, minshear]=Analyze_Data2_3D(w,dt,e0,stress,strain,numtimesteps)

%determine how many complete periods there are.
%first half period is the initial strain, should not be counted

[t_s,max_strain]=max(strain(100:end));
[~,min_strain]=min(strain(100:end));

[max_stress,t_max]=max(stress(100:end));
[min_stress,t_min]=min(stress(100:end));

    
Delta=(t_max+1-max_strain)*dt;
% Delta=mod(Delta,2*pi);

% s=mean(max_stress-min_stress);

maxshear=mean(max_stress);
minshear=mean(min_stress);

G1=max_stress/t_s*cos(Delta);
G2=max_stress/t_s*sin(Delta);


   