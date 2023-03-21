clc
clear all
% close all;

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
% Tstart=0.0;
% Rs=0.0;
% Vdcref=270;
% alpha=1e10;         %% 0.5
% Cdc=150e-6; 
% R=90*alpha/(1+alpha);
% RCdc=0.070;
% P=Vdcref^2/R;
% Vse=57.5;
% Vse=Vse+P/3/Vse*Rs*1;
% Vsm=Vse*sqrt(3);
% 
% fline=400; w=2*pi*fline;
% fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
% 
% %% low pass filter
% k=1/sqrt(2);
% wf=k*fline*2*pi;
% tf_LP = tf([0 wf],[1 wf]);
% tf_LP_d = c2d(tf_LP,1/20e3)
% 
% %% pll PI
% DEF_pll=50;							%control loop for PLL (Hz)
% DEF_pll_damp=0.707*1;%0.707;					%damping factor for the PLL controller
% DEF_Vin=Vse;						%input phase to neutral rms voltage
% DEF_Tsw=Tsw;						%switching period
% FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
% FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
% 
% tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
% PI_pll=d2c(tf_pll_z)

s= tf([1 0],[0 1]);
PI_pll = 3.2 + 0.05/s;
%% integrator

I_pll = tf([0 1],[1 0]);
%% PLL loop gain
Tpll = 216.8*sqrt(3)*PI_pll*I_pll;
% Tpll_1 = 57.1*tf_LP*PI_pll*I_pll;
figure (1)
bode(Tpll)
hold on
grid on

