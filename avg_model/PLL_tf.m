clc
clear all
close all
%%
DEF_pll=10;%200;							%control loop for PLL (Hz)
DEF_pll_damp=0.4;%0.707;					%damping factor for the PLL controller
DEF_Vin=57.5;						%input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;

tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z)
I_pll = tf([0 1],[1 0]);
Tpll = 57.5*sqrt(2)*tf_pll*I_pll;
% kp_pll = 0.05;
% ki_pll = 3.2*10;
% tf_pll_d = tf([kp_pll ki_pll],[1 0])


figure(1)
bode(tf_pll,Tpll)
hold on
% bode(tf_pll_d)
% bode(tf_pll_d*2)