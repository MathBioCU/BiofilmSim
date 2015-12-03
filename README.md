# BiofilmSim
Simulation software of biofilms based on our extension to the Immersed Boundary Method.
To run simulations, use RunSimulation.m. Requires access to data file for initial bacteria coordinates. 

Note: the Analyze_Data2D does not return the right results. Instead, use the following matlab commands:

[t,m]=max(fStrain); %or min
[t2,m2]=max(v0*visc0/charLength*vShear-eShear); %or min

delta=(t-t2)*w;
G1=m2/m*cos(delta);
G2=m2/m*sin(delta);

