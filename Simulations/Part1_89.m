clc
close all
clear

G = tf([-0.06993 0 24.5], [1 1.713+6.39*3.496 0.7577+58.78*3.496 0 -132.4]);
Gcontrol = tf([304.4 311.8], [1 12.48])
ControlledSystem = G * Gcontrol;
G_final = feedback(ControlledSystem,1);
detail3 = stepinfo(G_final)
step(G_final)
%K >0
figure
rlocus(ControlledSystem)
title("Root Locus (K>0)")
%K <0
figure
rlocus(-ControlledSystem)
title("Root Locus (K<0)")