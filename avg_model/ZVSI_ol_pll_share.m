clc
clear all;
close all

Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1e0 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=0e-6;
Vdcref=270;
alpha=1e10;         %% 0.5
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=1000e-6; 
RL=110e-3;
f=400; omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
s=tf([1 0],[0 1]);


Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);

Gsvm = [Gsvmdd 0*Gsvmdd;-0*Gsvmdd Gsvmdd];

fprintf('initialization is done!\n')


I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Gdel = Gsvm;


f_pll = [50 100 200];
Dd_ref = 0.3734;
Dq_ref = 0.1027;
for i=1:1:1
%% PLL
DEF_pll=f_pll(2);							%control loop for PLL (Hz)
DEF_pll_damp=0.707;                         %damping factor for the PLL controller
DEF_Vin=57.5;                               %input phase to neutral rms voltage
DEF_Tsw=0.00005;                            %switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z);

Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd=Vsdq(1);
Vsq=Vsdq(2);
Dd=Dd_ref;
Dq=Dq_ref;
Id=-11.03;
Iq=0;
fprintf('initialization is done!\n')

%% Calculation based on model
% s=tf([1 0],[0 1]);
% Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
% omega=w;
Vdc=270;
Zo_cal=[GL -1*omega*L; 1*omega*L GL];
Yin_cal=1/Zo_cal;


den1=(GL)*(GL)-(-1*omega*L)*(1*omega*L);
num11=(Vdc)*(GL);
num12=(Vdc)*(1*omega*L);
num22=num11;
num21=-num12;
Gid_cal=[-num11/den1 -num12/den1;-num21/den1 -num22/den1];

E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[0 Iq*Gpll;0 -Id*Gpll];
Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];


Zin_pll_cal=1/(1/Zo_cal+Gid_cal*Gdel*Gdpll);
Yin_cal=1/Zo_cal;
Ypll=Gid_cal*Gdel*Gdpll;
Yin_pll_cal=Yin_cal+Ypll;
%% calculation and simulation comparison
figure(1)
bode(Zo_cal,Bode_O)
legend('Zin\_cal')
Bode_Darklines(3)

figure(3)
bode(Gid_cal,Bode_O)
legend('Gid\_cal')
Bode_Darklines(3)

figure(5)
bode(Zin_pll_cal,Bode_O)
legend('Zin\_pll\_cal')
Bode_Darklines(3)

Zin_pll_1_cal=Zin_pll_cal;


%% plot open loop and open loop with pll
figHandle=figure(11)
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zo_cal,Zin_pll_cal,Bode_O)
Bode_Darklines(3)

end