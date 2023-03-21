% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFE parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd ..
% cd ..
% 
% str=strcat(pwd,'\VICTO_simulation_files');
% addpath(str)
% 
% str=strcat(pwd,'\My Libraries\Virtual Impedance Analyzer');
% addpath(str)
% cd Models\sw_model_SC

clear function
Vse=57.5;
Rss=0.5;
Vsm=Vse*sqrt(2)-8.2/sqrt(2)*Rss;
%% boost inductor parameters
Lboost_a    = 474e-6;           % inductor No. 2
RLboost_a   = 0.080;
Lboost_b    = 471e-6;           % inductor No. 1
RLboost_b   = 0.084;
Lboost_c    = 463e-6;           % inductor No. 3
RLboost_c   = 0.110;
%% dc link cap parameters
Cdc_afe     = 98.8e-6; 
RCdc_afe    = 0.049;
LCdc_afe    = 0.65e-6;
%% load resistor
Rdc         = 96;
%% dead time and PWM parameters
DeadTime =0.0e-6;
PWM_Cycle=4000;
%% system parameters
f=400; w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
m=0.4;
%parameters for IGBT 6MBP30RH060-50
Vfs=1.0;
Vfd=0.75;
Rsd=0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VSI parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output filter inductor parameters
La  = 1000e-6;%962e-6;       % inductor No. 1
RLa = 0.100;        
Lb  = 1000e-6;%985e-6;       % inductor No. 2
RLb = 0.100;
Lc  = 1000e-6;%972e-6;       % inductor No. 3
RLc = 0.100;
%% output filter cap parameters
Ca  = 32e-6;        % cap No. 1
LCa = 0.61e-6;
RCa = 0.051;
Cb  = 31.5e-6;      % cap No. 2
LCb = 0.61e-6;
RCb = 0.054;
Cc  = 31.8e-6;      % cap No. 3
LCc = 0.63e-6;
RCc = 0.052;
%% dc link cap parameters
Cdc_vsi     = 147.3e-6;
LCdc_vsi    = 0.57e-6;
RCdc_vsi    = 0.046;
%% load resistor
R=16;
%% input dc voltage
Vdc=270/1;
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);

Pref_VSI    = -0e3;
Qref_VSI    = -1e3;

Kpi_VSI		= 6.2832/Vdc;
Kii_VSI		= 628.3185/Vdc;
Kppq_VSI	= 70*4e-5;
Kipq_VSI	= 70;
%% VSI PLL parameters
DEF_pll=100;							%control loop for PLL (Hz)
DEF_pll_damp=0.707;                     %damping factor for the PLL controller
DEF_Vin=57.5;                           %input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z)
    
Kp_PLL_VSI  = 8.921;
Ki_PLL_VSI  = 3964;

Kpi_AFE		= 0.0116*270;
Kii_AFE		= 46.5421*270;
Kpv_AFE		= 0.0462;
Kiv_AFE		= 4.5815;

Kp_PLL_AFE  = 0.2;
Ki_PLL_AFE  = 1;

fprintf('\ndone\n')