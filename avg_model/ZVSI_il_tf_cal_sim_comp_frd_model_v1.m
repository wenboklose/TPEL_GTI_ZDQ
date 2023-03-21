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
Rs=0*1e-0;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=105e-6; 
R=96;
RCdc=0.049;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=1000e-6;
RL=110e-3;
I = [1 0; 0 1];
f=400;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2);
Id_ref_a=[-6.03 -11.03 -20.03];%-7.6907;%-11.03;%7.6907;
Iq_ref_a=[-6 0 6];
Vdc=270;
s=tf([1 0],[0 1]);
Ns= [0.0 0.0045 0.009];     %% modified coefficient for pll
f_pll = [50 100 200];
ki_i = [1 10 20];
for i=1:1:1
%     N=Ns(i);
    %% PLL
    Id_ref=1*Id_ref_a(2);
    Iq_ref=1*Iq_ref_a(2);
    DEF_pll=f_pll(2);					%control loop for PLL (Hz)
    DEF_pll_damp=0.707;					%damping factor for the PLL controller
    DEF_Vin=57.5;						%input phase to neutral rms voltage
    DEF_Tsw=0.00005;					%switching period
    FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
    FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
    tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
    tf_pll=d2c(tf_pll_z);
%     Kp1 = 1;%0.1; %No.1 PLL Kp (inverter)
%     Ki1 = 2; %2; %No.1 PLL Ki 3.2 is unstable
%     tf_pll = Kp1+Ki1/s;
    
    
    fv=100; fi=1000;
    Lboost_con=1000e-6;
    RLboost_con=110e-3;
    Cdc_con=100e-6;
    Rdc_con=90;
    kpv=2*pi*fv*Cdc_con;%4.4396*0.06;%2*pi*fv*Cdc_con; %0.3770;%
    kiv=2*pi*fv/Rdc_con;%4.4396;%2*pi*fv/Rdc_con; %27.9253;%
    kpi=2*pi*fi*Lboost_con/Vdcref; %0.0233;% ki_i(2)*0.00068;%
    kii=ki_i(2)*2*pi*fi*RLboost_con/Vdcref; %4.6542;%ki_i(2);%
% 
%     Vsdq=[sqrt(3/2)*Vsm; 0];
%     Vsd=Vsdq(1);
%     Vsq=Vsdq(2);
% 
%     Id=7.6907;
%     Iq=0.00;

    %% signal conditioning filter
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
%     Gsvm = [Gsvm(1,1) 0*Gsvm(1,2);-0*Gsvm(2,1) Gsvm(2,2)];
    Gsvm = [Gsvmdd 0*Gsvm(1,2);-0*Gsvm(2,1) Gsvmdd];
    
    fprintf('initialization is done!\n')

    % sys = linearize('AFE_avg_il_pll',0.5);
    I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
    Gdel = Gsvm;
    Ki = tf_filter_dq;
    Kv = tf_filter_dq;%Hv;
    open('VSI_GT_avg_il_pll_sv.slx')
    sys = linearize('VSI_GT_avg_il_pll_sv',0.5);
    %% Zin model linearization
    H=ss(sys);
    Vavg1=[H(1,1) H(1,2);
        H(2,1) H(2,2);];
    Ilavg1=[H(3,1) H(3,2);
        H(4,1) H(4,2);];
    Zin_il_pll_avg_sim=Vavg1/Ilavg1;

    Dd = Dd_rec(length(Dd_rec)-10)
    Dq = Dq_rec(length(Dq_rec)-10)
    Id = Id_rec(length(Id_rec)-10);
    Iq = 1*Iq_rec(length(Iq_rec)-10)
    Vsd = Vsd_rec(length(Vsd_rec)-10);
    Vsq = Vsq_rec(length(Vsq_rec)-10);
    Vdc = Vdc_rec(length(Vdc_rec)-10);


    %% Calculation based on model of frd
    n= 40;
    w = logspace(-4,4,n);
    w = 2*pi*w;

    E0=Vsd;
    Gpll=tf_pll/(s+E0*tf_pll);

    Gipll=[tf([0 0],[0 1]) Iq*Gpll;tf([0 0],[0 1]) -Id*Gpll];

    Gdel_frd = [frd(freqresp(Gdel(1,1),w),w) frd(freqresp(Gdel(1,2),w),w);...
        frd(freqresp(Gdel(2,1),w),w) frd(freqresp(Gdel(2,2),w),w)];

    Ki_frd = [frd(freqresp(Ki(1,1),w),w) frd(freqresp(Ki(1,2),w),w);...
        frd(freqresp(Ki(2,1),w),w) frd(freqresp(Ki(2,2),w),w)];

    Kv_frd = [frd(freqresp(Kv(1,1),w),w) frd(freqresp(Kv(1,2),w),w);...
        frd(freqresp(Kv(2,1),w),w) frd(freqresp(Kv(2,2),w),w)];

    Gipll_frd = [frd(freqresp(Gipll(1,1),w),w) frd(freqresp(Gipll(1,2),w),w);...
        frd(freqresp(Gipll(2,1),w),w) frd(freqresp(Gipll(2,2),w),w)];

    Gdpll=[tf([0 0],[0 1]) -Dq*Gpll;tf([0 0],[0 1]) +Dd*Gpll];
    Gdpll_frd = [frd(freqresp(Gdpll(1,1),w),w) frd(freqresp(Gdpll(1,2),w),w);...
        frd(freqresp(Gdpll(2,1),w),w) frd(freqresp(Gdpll(2,2),w),w)];
    GL=1*L*s+1*RL;
    
    Zo_cal=[GL -1*omega*L; 1*omega*L GL];
    Yin_cal=1/Zo_cal;
    Yin_cal_frd = [frd(freqresp(Yin_cal(1,1),w),w) frd(freqresp(Yin_cal(1,2),w),w);...
        frd(freqresp(Yin_cal(2,1),w),w) frd(freqresp(Yin_cal(2,2),w),w)];

    Dt=[Dd Dq;0 0];
    It=[Id Iq;0 0];
%     Gve_cal=Zdc*Dt*Yin_cal;
%     Gve_cal_frd = [frd(freqresp(Gve_cal(1,1),w),w) frd(freqresp(Gve_cal(1,2),w),w);...
%         frd(freqresp(Gve_cal(2,1),w),w) frd(freqresp(Gve_cal(2,2),w),w)];
    den1=(GL)*(GL)-(-1*omega*L)*(1*omega*L);
    num11=(Vdc)*(GL);
    num12=(Vdc)*(1*omega*L);
    num22=num11;
    num21=-num12;
    Gid_cal=[-num11/den1 -num12/den1;-num21/den1 -num22/den1];
    Gid_cal_frd = [frd(freqresp(Gid_cal(1,1),w),w) frd(freqresp(Gid_cal(1,2),w),w);...
        frd(freqresp(Gid_cal(2,1),w),w) frd(freqresp(Gid_cal(2,2),w),w)];
%     Gvd_cal=Zdc*(Dt*Gid_cal+1*It);
%     Gvd_cal_frd = [frd(freqresp(Gvd_cal(1,1),w),w) frd(freqresp(Gvd_cal(1,2),w),w);...
%         frd(freqresp(Gvd_cal(2,1),w),w) frd(freqresp(Gvd_cal(2,2),w),w)];

    % Ypll_frd=Gid_cal_frd*Gdpll_frd;
    % Yin_pll_cal_frd=Yin_cal_frd+Ypll_frd;

    %% Calculatoin with current control loop
    I_frd = [frd(freqresp(I(1,1),w),w) frd(freqresp(I(1,2),w),w);...
        frd(freqresp(I(2,1),w),w) frd(freqresp(I(2,2),w),w)];
    Gdei = [tf([0 0],[0 1]) tf([0 1/Vdcref*Lboost_con*omega],[0 1]);...
        -tf([0 1/Vdcref*Lboost_con*omega],[0 1]) tf([0 0],[0 1])];
%     Gdei = [tf([0 0],[0 1]) tf([0 0],[0 1]);...
%         tf([0 0],[0 1]) tf([0 0],[0 1])];
    Gdei_frd = [frd(freqresp(Gdei(1,1),w),w) frd(freqresp(Gdei(1,2),w),w);...
        frd(freqresp(Gdei(2,1),w),w) frd(freqresp(Gdei(2,2),w),w)];
    Gci = -[kpi+kii/s tf([0 0],[0 1]); tf([0 0],[0 1]) kpi+kii/s];
    Gci_frd = [frd(freqresp(Gci(1,1),w),w) frd(freqresp(Gci(1,2),w),w);...
        frd(freqresp(Gci(2,1),w),w) frd(freqresp(Gci(2,2),w),w)];

    Ypll_il_frd = (-Gid_cal_frd*Gdel_frd*((-Gdei_frd + Gci_frd)*Gipll_frd*Kv_frd-Gdpll_frd*Kv_frd));
    T_frd = (I_frd+Gid_cal_frd*Gdel_frd*(-Gdei_frd + Gci_frd)*Ki_frd);
    Yin_il_pll_avg_cal_frd = T_frd\(Ypll_il_frd+Yin_cal_frd);
    Zin_il_pll_avg_cal_frd = (Ypll_il_frd+Yin_cal_frd)\T_frd;
    Ti_frd = Gid_cal_frd*Gdel_frd*(-Gdei_frd + Gci_frd)*Ki_frd;

    figHandle=figure(5)
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Zin_il_pll_avg_sim,Zin_il_pll_avg_cal_frd,Bode_O)
%     bode(Zin_il_pll_avg_cal_frd,Bode_O)
%     legend('Zin\_il\_pll','Y\_PLL\_o','Yin\_ol\_PLL\_cal')
    hold on
    Bode_Darklines(3)
    
    figure(51)
    bode(Zin_il_pll_avg_sim,Bode_O)
%     legend('Zin\_il\_pll','Y\_PLL\_o','Yin\_ol\_PLL\_cal')
    hold on
    Bode_Darklines(3)
%     
    figHandle=figure(6)
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    pzmap(Zin_il_pll_avg_sim)
    hold on

    %% Zdd and Zqq simplified equations:
    Zdds = (L*s^2+(RL+Vdc*kpi)*s+Vdc*kii)/s;
    Zdds = Vdc*kii/s;
    Zqqs = 1/(1/Zdds-1/Zdds*Vsd*Gpll+Id*Gpll);
    kppll = 8.921; kipll = 3964;
    c = Id*Vdc*kii*kppll; d = Id*Vdc*kii*kipll;
    Root1 = (-d/2+sqrt(d^2/4+c^3/27))^(1/3)+(-d/2-sqrt(d^2/4+c^3/27))^(1/3);
    Root1_p1 = (-d/2+sqrt(d^2/4+c^3/27))^(1/3);
    Root1_p2 = (-d/2-sqrt(d^2/4+c^3/27))^(1/3);
    figHandle=figure(7)
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Zin_il_pll_avg_cal_frd(1,1),Zin_il_pll_avg_cal_frd(1,2),Zin_il_pll_avg_cal_frd(2,1),Zin_il_pll_avg_cal_frd(2,2),Bode_O)
    hold on
    Bode_Darklines(3)
    figHandle=figure(71)
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Zin_il_pll_avg_cal_frd(1,1),Zdds,Bode_O)
    hold on
    Bode_Darklines(3)
    figHandle=figure(72)
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Zin_il_pll_avg_cal_frd(2,2),Zqqs,Bode_O)
    hold on
    Bode_Darklines(3)
end
% bode(Zo_cal,Bode_O)
% Bode_Darklines(3)
% cd ..
% cd exp_measurement
% load Z_PLL_50_Id_11_3.mat
% Z_PLL_50_11_03_exp = Z.ZDQRAW/3;
% Y_PLL_50_11_03_exp = 1/Z_PLL_50_11_03_exp;
% load Z_PLL_100_Id_11_3.mat
% Z_PLL_100_11_03_exp = Z.ZDQRAW/3;
% Y_PLL_100_11_03_exp = 1/Z_PLL_100_11_03_exp;
% load Z_PLL_200_Id_11_3.mat
% Z_PLL_200_11_03_exp = Z.ZDQRAW/3;
% Y_PLL_200_11_03_exp = 1/Z_PLL_200_11_03_exp;
% load Z_PLL_100_Id_11_3.mat
% Z_PLL_100_6_03_exp = Z.ZDQRAW/3;
% Y_PLL_100_6_03_exp = 1/Z_PLL_100_6_03_exp;
% load Z_PLL_100_Id_14_3.mat
% Z_PLL_100_14_03_exp = Z.ZDQRAW/3;
% Y_PLL_100_14_03_exp = 1/Z_PLL_100_14_03_exp;
% cd ..
% cd avg_model
% 
% figHandle=figure(6)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Zin_il_pll_avg_cal_frd,Z_PLL_200_11_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% figHandle=figure(7)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Yin_il_pll_avg_cal_frd,Y_PLL_200_11_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% 
% 
% figHandle=figure(8)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Zin_il_pll_avg_cal_frd,Z_PLL_100_11_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% figHandle=figure(9)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Yin_il_pll_avg_cal_frd,Y_PLL_100_11_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% 
% 
% figHandle=figure(10)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Zin_il_pll_avg_cal_frd,Z_PLL_50_11_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% figHandle=figure(11)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Yin_il_pll_avg_cal_frd,Y_PLL_50_11_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% 
% figHandle=figure(12)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Zin_il_pll_avg_cal_frd,Z_PLL_100_6_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% figHandle=figure(13)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Yin_il_pll_avg_cal_frd,Y_PLL_100_6_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% 
% figHandle=figure(14)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Zin_il_pll_avg_cal_frd,Z_PLL_100_14_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
% figHandle=figure(15)
% set(figHandle, 'Position', [50, 10, 1000, 800]);
% bode(Yin_il_pll_avg_cal_frd,Y_PLL_100_14_03_exp,Bode_O)
% hold on
% Bode_Darklines(3)
