clc
clear all;
close all

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
    Id_ref=-11.03;
    Iq_ref=0;
    DEF_pll=100;	
    DEF_pll_damp=0.707;
    DEF_Vin=57.5;
    DEF_Tsw=0.00005;
    FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;	
    FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
    tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
    tf_pll=d2c(tf_pll_z);
    
    
    fv=100; fi=1000;
    Lboost_con=1000e-6;
    RLboost_con=110e-3;
    Cdc_con=100e-6;
    Rdc_con=90;
    kpv=0.0628;
    kiv=6.9813;
    kpi=0.0233;
    kii=25.5982;

    %% signal conditioning filter
    C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
    tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
    tf_filter_dq = [tf_filter 0*tf_filter; 0*tf_filter tf_filter];

    Tdelay = 1.5/fsw;
    Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);
    Gsvm = [Gsvmdd 0*Gsvmdd;-0*Gsvmdd Gsvmdd];
    
    Dd = 0.3734
    Dq = 0.1027
    Id = -11.0300
    Iq = 0
    Vsd = 99.5929
    Vsq = 0
    Vdc = 270


    %% Calculation based on model of frd
    n= 40;
    w = logspace(0,4,n);
    w = 2*pi*w;

    E0=Vsd;
    Gpll=tf_pll/(s+E0*tf_pll);

    Gipll=[tf([0 0],[0 1]) Iq*Gpll;tf([0 0],[0 1]) -Id*Gpll];
    Gdel = Gsvm;
    Ki = tf_filter_dq;
    Kv = tf_filter_dq;
    I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
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
    den1=(GL)*(GL)-(-1*omega*L)*(1*omega*L);
    num11=(Vdc)*(GL);
    num12=(Vdc)*(1*omega*L);
    num22=num11;
    num21=-num12;
    Gid_cal=[-num11/den1 -num12/den1;-num21/den1 -num22/den1];
    Gid_cal_frd = [frd(freqresp(Gid_cal(1,1),w),w) frd(freqresp(Gid_cal(1,2),w),w);...
        frd(freqresp(Gid_cal(2,1),w),w) frd(freqresp(Gid_cal(2,2),w),w)];
    %% Calculatoin with current control loop
    I_frd = [frd(freqresp(I(1,1),w),w) frd(freqresp(I(1,2),w),w);...
        frd(freqresp(I(2,1),w),w) frd(freqresp(I(2,2),w),w)];
    Gdei = [tf([0 0],[0 1]) tf([0 1/Vdcref*Lboost_con*omega],[0 1]);...
        -tf([0 1/Vdcref*Lboost_con*omega],[0 1]) tf([0 0],[0 1])];
    Gdei_frd = 0*[frd(freqresp(Gdei(1,1),w),w) frd(freqresp(Gdei(1,2),w),w);...
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
    bode(Zin_il_pll_avg_cal_frd)
    
    
end
