load('../../DynamicModuliData.mat');

figure(1)
P1=loglog(AngFreq,G1_86,'-b');
hold on
P2=loglog(AngFreq,G2_86,'-g');

P3=loglog(AngFreq,G1_86+G1_86e,'--.k');
P4=loglog(AngFreq,G1_86-G1_86e,'--.k');
P5=loglog(AngFreq,G2_86+G2_86e,'--.k');
P6=loglog(AngFreq,G2_86-G2_86e,'--.k');
