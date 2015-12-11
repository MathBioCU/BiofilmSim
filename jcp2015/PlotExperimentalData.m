load('../../DynamicModuliData.mat');

figure(1)
loglog(AngFreq,G1_86,'-');
hold on
loglog(AngFreq,G2_86,'-.');

loglog(AngFreq,G1_86+G1_86e,'--.k');
loglog(AngFreq,G1_86-G1_86e,'--.k');
loglog(AngFreq,G2_86+G2_86e,'--.k');
loglog(AngFreq,G2_86-G2_86e,'--.k');
