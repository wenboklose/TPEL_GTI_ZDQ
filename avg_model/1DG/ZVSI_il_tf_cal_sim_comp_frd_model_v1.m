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
Bode_O.XLim={[1e-3 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
Tstart=0.015;
Rs=1e-1;
Vdcref=600;
alpha=1e10;         %% 0.5
Cdc=105e-6; 
R=8.87;%8.964;%8.365;
RCdc=0.049;
P=Vdcref^2/R;
Vse=216.8/sqrt(2);
% Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=1000e-6;
RL=110e-3;
I = [1 0; 0 1];
f=60;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2);
% Id_ref=-11.03;      % 8 A Id current sinking by AFE for Dong Dong's system; -11.03 A Id current given by VSI
id_ref=-190;
Iq_ref=0;
Vdc=600;
s=tf([1 0],[0 1]);
% Ns= [0.0 0.0045 0.009];     %% modified coefficient for pll
% f_pll = [10 20 40];
% ki_i = [20 30 40];
Ki_pll = [3.2 32];
for i=1:1:1
    Id_ref = id_ref;
    Kp1 = 1.5;%3;% 1.5 is stable case, 3 is unstable case
    Ki1 = Ki_pll(1); %this value is fixed for 3.2
    tf_pll = Kp1+Ki1/s;
    
    fv=100; fi=1000;
    Lboost_con=1000e-6;
    RLboost_con=110e-3;
    Cdc_con=100e-6;
    Rdc_con=90;
    kpv=2*pi*fv*Cdc_con;%4.4396*0.06;%2*pi*fv*Cdc_con; %0.3770;%
    kiv=2*pi*fv/Rdc_con;%4.4396;%2*pi*fv/Rdc_con; %27.9253;%
    kpi=2*pi*fi*Lboost_con/Vdcref; %0.0233;% ki_i(2)*0.00068;%
    kii=2*pi*fi*RLboost_con/Vdcref; %4.6542;%ki_i(2);%
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
    Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];

    fprintf('initialization is done!\n')

    % sys = linearize('AFE_avg_il_pll',0.5);
    I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
    Gdel = Gsvm;
    Ki = tf_filter_dq;
    Kv = tf_filter_dq;%Hv;
    open('VSI_GT_avg_il_pll.mdl')
    sys = linearize('VSI_GT_avg_il_pll',.3);
    %% Zin model linearization
    H=ss(sys);
    Vavg1=[H(1,1) H(1,2);
        H(2,1) H(2,2);];
    Ilavg1=[H(3,1) H(3,2);
        H(4,1) H(4,2);];
    Zin_il_pll_avg_sim=Vavg1/Ilavg1;

    Dd = Dd_rec(length(Dd_rec)-10);
    Dq = Dq_rec(length(Dq_rec)-10);
    Id = Id_rec(length(Id_rec)-10);
    Iq = 1*Iq_rec(length(Iq_rec)-10);
    Vsd = Vsd_rec(length(Vsd_rec)-10);
    Vsq = Vsq_rec(length(Vsq_rec)-10);
    Vdc = Vdc_rec(length(Vdc_rec)-10);


    %% Calculation based on model of frd
    n= 1e3;
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
    % Gdei = [tf([0 0],[0 1]) tf([0 0],[0 1]);...
    %     tf([0 0],[0 1]) tf([0 0],[0 1])];
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


    figure(5)
%     bode(Zin_il_pll_avg_cal_frd,Zin_il_pll_avg_sim,Bode_O)
    bode(Zin_il_pll_avg_sim,Bode_O)
    hold on
    Bode_Darklines(3)
    figure(6)
    bode(1/Zin_il_pll_avg_cal_frd,1/Zin_il_pll_avg_sim,Bode_O)
    hold on
    Bode_Darklines(3)

    
    figure(7)
%     bode(Zin_il_pll_avg_cal_frd,Zin_il_pll_avg_sim,Bode_O)
    bode(Zin_il_pll_avg_sim(2,2),Bode_O)
    hold on
    Bode_Darklines(3)
    
end
% bode(Zo_cal,Bode_O)
% Bode_Darklines(3)
