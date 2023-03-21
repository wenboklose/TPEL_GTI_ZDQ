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
%% on board signal condition filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); 0*tf_filter_dq(2,1) tf_filter_dq(2,2)];

Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);

%     Gsvm = [Gsvmdd*cos(omega*Tdelay) 0*Gsvmdd*sin(omega*Tdelay);-0*Gsvmdd*sin(omega*Tdelay) Gsvmdd*cos(omega*Tdelay)];
Gsvm = tf(JF_DQFromABC(Gsvmdd,omega));
Gsvm = [Gsvm(1,1) 0*Gsvm(1,2);-0*Gsvm(2,1) Gsvm(2,2)];
Gsvm = [Gsvmdd 0*Gsvm(1,2);-0*Gsvm(2,1) Gsvmdd];

fprintf('initialization is done!\n')

% sys = linearize('AFE_avg_il_pll',0.5);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Gdel = Gsvm;
Ki = tf_filter_dq;
Kv = tf_filter_dq;%Hv;

f_pll = [50 100 200];
Dd_ref = [0.3713 0.3734 0.3770];
Dq_ref = [0.0561 0.1027 0.1864];
for i=1:1:1
%% PLL
DEF_pll=f_pll(2);							%control loop for PLL (Hz)
DEF_pll_damp=0.707;					%damping factor for the PLL controller
DEF_Vin=57.5;						%input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z);

% fv=100; fi=1000;
% Lboost_con=500e-6;
% RLboost_con=100e-3;
% Cdc_con=150e-6;
% Rdc_con=90;
% kpv=2*pi*fv*Cdc_con; kiv=2*pi*fv/Rdc_con;
% kpi=2*pi*fi*Lboost_con/Vdcref; kii=2*pi*fi*RLboost_con/Vdcref;

Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd=Vsdq(1);
Vsq=Vsdq(2);
Dd=Dd_ref(2);%0.3803;%
Dq=Dq_ref(2);%0.0154;%
% Id=-11.03;
% Iq=0;
fprintf('initialization is done!\n')


% %% Gvd model linearization:
% model = 'AFE_avg_Gvd';
% 
% %% Create the linearization I/O as specified in AFE_avg_Gvd
% ios(3) = linio('AFE_avg_Gvd/vdc',1,'out');
% ios(2) = linio('AFE_avg_Gvd/Dq',1,'in');
% ios(1) = linio('AFE_avg_Gvd/Dd',1,'in');
% 
% %% Linearize the model
% Gvd_avg_sim = linearize(model,0.5,ios);
% 
% %% Gve model linearization:
% model = 'AFE_avg_Gve';
% 
% %% Create the linearization I/O as specified in AFE_avg_Gve
% ios(3) = linio('AFE_avg_Gve/vdc',1,'out');
% ios(2) = linio('AFE_avg_Gve/vpq',1,'in');
% ios(1) = linio('AFE_avg_Gve/vpd',1,'in');
% 
% %% Linearize the model
% Gve_avg_sim = linearize(model,0.5,ios);

%% Gid model linearization:
model = 'VSI_GT_avg_Gid';

%% Create the linearization I/O as specified in AFE_avg_Gid
ios(4) = linio('VSI_GT_avg_Gid/ilq',1,'out');
ios(3) = linio('VSI_GT_avg_Gid/ild',1,'out');
ios(2) = linio('VSI_GT_avg_Gid/Dq',1,'in');
ios(1) = linio('VSI_GT_avg_Gid/Dd',1,'in');

%% Linearize the model
Gid_avg_sim = linearize(model,0.5,ios);

%% Zin model linearization:
model = 'VSI_GT_avg_Zol';

sys = linearize('VSI_GT_avg_Zol',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_avg_sim=Vavg1/Ilavg1;

%% Zin_pll model linearization:
model = 'VSI_GT_avg_Zol_pll';

sys = linearize('VSI_GT_avg_Zol_pll',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zin_pll_avg_sim=Vavg1/Ilavg1;
Dd = Dd_rec(length(Dd_rec)-10)
    Dq = Dq_rec(length(Dq_rec)-10)
    Id = Id_rec(length(Id_rec)-10);
    Iq = 1*Iq_rec(length(Iq_rec)-10);
    Vsd = Vsd_rec(length(Vsd_rec)-10);
    Vsq = Vsq_rec(length(Vsq_rec)-10);
    Vdc = Vdc_rec(length(Vdc_rec)-10);
%% Calculation based on model
% s=tf([1 0],[0 1]);
% Zdc=R/(R*Cdc*s+1);
GL=1*L*s+RL;
% omega=w;
Vdc=270;
Zo_cal=[GL -1*omega*L; 1*omega*L GL];
Yin_cal=1/Zo_cal;
% Dt=[Dd Dq;0 0];
% It=[Id Iq;0 0];

den1=(GL)*(GL)-(-1*omega*L)*(1*omega*L);
num11=(Vdc)*(GL);
num12=(Vdc)*(1*omega*L);
num22=num11;
num21=-num12;
Gid_cal=[-num11/den1 -num12/den1;-num21/den1 -num22/den1];




% Vse=57.5;
% Vsm=Vse*sqrt(2);
% Vsdq=[sqrt(3/2)*Vsm; 0];
% Vsd=Vsdq(1)
% Vsq=Vsdq(2)
E0=Vsd;
Gpll=tf_pll/(s+E0*tf_pll);
Gipll=[0 Iq*Gpll;0 -Id*Gpll];
Gdpll=[0 -Dq*Gpll;0 +Dd*Gpll];

% Zin_ol_pll=1/(Gppll+Yin_ol*Gspll);
Zin_pll_cal=1/(1/Zo_cal+Gid_cal*Gdel*Gdpll*Kv);
Yin_cal=1/Zo_cal;
Ypll=Gid_cal*Gdel*Gdpll*Kv;
Yin_pll_cal=Yin_cal+Ypll;
%% calculation and simulation comparison
figure(1)
bode(Zin_avg_sim,Zo_cal,Bode_O)
legend('Zin\_avg\_sim','Zin\_cal')
Bode_Darklines(3)

% figure(2)
% bode(Gve_avg_sim,Gve_cal(1,:),Bode_O)
% legend('Gve\_avg\_sim','Gve\_cal')
% Bode_Darklines(3)

figure(3)
bode(Gid_avg_sim,Gid_cal,Bode_O)
legend('Gid\_avg\_sim','Gid\_cal')
Bode_Darklines(3)

% figure(4)
% bode(Gvd_avg_sim,Gvd_cal(1,:),Bode_O)
% legend('Gvd\_avg\_sim','Gvd\_cal')
% Bode_Darklines(3)

figure(5)
bode(Zin_pll_avg_sim,Zin_pll_cal,Bode_O)
legend('Zin\_pll\_avg\_sim','Zin\_pll\_cal')
Bode_Darklines(3)
Zin_pll_1_avg_sim=Zin_pll_avg_sim;
Zin_pll_1_cal=Zin_pll_cal;
% save('ZAFE_in_pll_800.mat','Zin_pll_1_avg_sim','Zin_pll_1_cal');

figure(6)
bode(Yin_cal,Ypll,Yin_pll_cal,Bode_O)
legend('Yin','Y\_PLL\_o','Yin\_ol\_PLL\_cal')
Bode_Darklines(3)

figure(8)
bode(Yin_cal(2,2),Ypll(2,2),Yin_pll_cal(2,2),Bode_O)
legend('Yin\_qq','Y\_PLL\_o\_qq','Yin\_ol\_PLL\_cal\_qq')
Bode_Darklines(3)

figure(7)
bode(Gid_cal,Ypll,Yin_cal,Bode_O)
Bode_Darklines(3)

figure(9)
pzmap(Yin_cal,Ypll,Yin_cal+Ypll)

figure(10)
pzmap(Yin_cal(2,2),Ypll(2,2),Yin_cal(2,2)+Ypll(2,2))
figure(101)
pzmap(Yin_cal(2,2),Ypll(2,2))
hold on

figHandle=figure(102)
set(figHandle, 'Position', [50, 10, 1000, 800]);
pzmap(Zin_pll_cal(1,1),Zin_pll_cal(2,2))
hold on
figHandle=figure(103)
set(figHandle, 'Position', [50, 10, 1000, 800]);
pzmap(Zin_pll_cal(2,1),Zin_pll_cal(1,2))
hold on
figHandle=figure(104)
set(figHandle, 'Position', [50, 10, 1000, 800]);
pzmap(Zin_pll_cal)
hold on

%% plot open loop and open loop with pll
figHandle=figure(11)
set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Zin_avg_sim,Zo_cal,Zin_pll_avg_sim,Zin_pll_cal,Bode_O)
bode(Zo_cal,Zin_pll_cal,Bode_O)
Bode_Darklines(3)

figHandle=figure(12)
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Zin_pll_cal,Bode_O)
hold on
Bode_Darklines(3)
end