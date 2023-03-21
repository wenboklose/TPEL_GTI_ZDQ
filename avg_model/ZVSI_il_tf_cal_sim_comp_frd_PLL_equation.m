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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth=3; fontsize=18;
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='m';
linestyle5='k'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3];

s=tf([1 0],[0 1]);
% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
%% RLC passive load
RS = 2;
RL = 10;
CL = 400e-6;
LL = 17.4e-3;
RLC = 1/(1/RS+1/RL+1/(LL*s)+CL*s);
Z_RLC_dq = JF_DQFromABC(RLC,60*2*pi);
%% converter parameters
Tstart=0.015;
Rs=1e-1;
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
f=60;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2);
Id_ref=-10.04;%7.6907;
Iq_ref=0;
Vdc=270;
s=tf([1 0],[0 1]);
Ns= [0.000 0.005 0.015];     %% modified coefficient for pll  // for Kp1=1; Ki1=2
% Ns= [0.000 0.005 0.0065];     %% modified coefficient for pll // for DEF_PLL
f_pll = [20 30 40];
ki_i = [20 30 40];
for i=1:1:1
%     N=Ns(i);
    %% PLL
    DEF_pll=f_pll(3);							%control loop for PLL (Hz)
    DEF_pll_damp=0.707;                       %damping factor for the PLL controller
    DEF_Vin=57.5;                             %input phase to neutral rms voltage
    DEF_Tsw=0.00005;                          %switching period
    FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
    FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
    tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
    tf_pll=d2c(tf_pll_z);
%     Kp1 = 1;%0.1; %No.1 PLL Kp (inverter)
%     Ki1 = 2; %2; %No.1 PLL Ki 3.2 is unstable
%     tf_pll = Kp1+Ki1/s;
    
    
    fv=100; fi=2000;
    Lboost_con=1000e-6;
    RLboost_con=110e-3;
    Cdc_con=100e-6;
    Rdc_con=90;
    kpv=2*pi*fv*Cdc_con;%4.4396*0.06;%2*pi*fv*Cdc_con; %0.3770;%
    kiv=2*pi*fv/Rdc_con;%4.4396;%2*pi*fv/Rdc_con; %27.9253;%
    kpi=2*pi*fi*Lboost_con/Vdcref; %0.0233;% ki_i(2)*0.00068;%
    kii=2*pi*fi*RLboost_con/Vdcref; %4.6542;%ki_i(2);%

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
    sys = linearize('VSI_GT_avg_il_pll',0.2);
    %% Zin model linearization
    H=ss(sys);
    Vavg1=[H(1,1) H(1,2);
        H(2,1) H(2,2);];
    Ilavg1=[H(3,1) H(3,2);
        H(4,1) H(4,2);];
    Zo_il_pll_avg_sim=Vavg1/Ilavg1;

    Dd = Dd_rec(length(Dd_rec)-10);
    Dq = Dq_rec(length(Dq_rec)-10);
    Id = Id_rec(length(Id_rec)-10);
    Iq = 1*Iq_rec(length(Iq_rec)-10);
    Vsd = Vsd_rec(length(Vsd_rec)-10);
    Vsq = Vsq_rec(length(Vsq_rec)-10);
    Vdc = Vdc_rec(length(Vdc_rec)-10);

    %% Calculation based on model of frd
    n= 2e3;
    w = logspace(-4,4,n);
    w = 2*pi*w;

    E0=Vsd;
    Gpll=tf_pll/(s+E0*tf_pll);

    Gipll=[tf([0 0],[0 1]) Iq*Gpll;tf([0 0],[0 1]) -Id*Gpll];
%     Gipll=[tf([0 0],[0 1]) Iq*((1-N*tf_pll*E0)*tf_pll/(s+E0*tf_pll)+N*tf_pll);tf([0 0],[0 1]) -Id*((1-N*tf_pll*E0)*tf_pll/(s+E0*tf_pll)+N*tf_pll)];

    Gdel_frd = [frd(freqresp(Gdel(1,1),w),w) frd(freqresp(Gdel(1,2),w),w);...
        frd(freqresp(Gdel(2,1),w),w) frd(freqresp(Gdel(2,2),w),w)];

    Ki_frd = [frd(freqresp(Ki(1,1),w),w) frd(freqresp(Ki(1,2),w),w);...
        frd(freqresp(Ki(2,1),w),w) frd(freqresp(Ki(2,2),w),w)];

    Kv_frd = [frd(freqresp(Kv(1,1),w),w) frd(freqresp(Kv(1,2),w),w);...
        frd(freqresp(Kv(2,1),w),w) frd(freqresp(Kv(2,2),w),w)];

    Gipll_frd = [frd(freqresp(Gipll(1,1),w),w) frd(freqresp(Gipll(1,2),w),w);...
        frd(freqresp(Gipll(2,1),w),w) frd(freqresp(Gipll(2,2),w),w)];

    Gdpll=[tf([0 0],[0 1]) -Dq*Gpll;tf([0 0],[0 1]) +Dd*Gpll];
%     Gdpll=[tf([0 0],[0 1]) -Dq*((1-N*tf_pll*E0)*tf_pll/(s+E0*tf_pll)+N*tf_pll);tf([0 0],[0 1]) +Dd*((1-N*tf_pll*E0)*tf_pll/(s+E0*tf_pll)+N*tf_pll)];
    Gdpll_frd = [frd(freqresp(Gdpll(1,1),w),w) frd(freqresp(Gdpll(1,2),w),w);...
        frd(freqresp(Gdpll(2,1),w),w) frd(freqresp(Gdpll(2,2),w),w)];
    GL=1*L*s+1*RL;
    
    Zo_cal=[GL -1*omega*L; 1*omega*L GL];
    Yin_cal=1/Zo_cal;
    Yin_cal_frd = [frd(freqresp(Yin_cal(1,1),w),w) frd(freqresp(Yin_cal(1,2),w),w);...
        frd(freqresp(Yin_cal(2,1),w),w) frd(freqresp(Yin_cal(2,2),w),w)];

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
    Ypll_il = (-Gid_cal*Gdel*((-Gdei + Gci)*Gipll*Kv-Gdpll*Kv));
    T_frd = (I_frd+Gid_cal_frd*Gdel_frd*(-Gdei_frd + Gci_frd)*Ki_frd);
    T = (I+Gid_cal*Gdel*(-Gdei + Gci)*Ki);
    Yin_il_pll_avg_cal_frd = T_frd\(Ypll_il_frd+Yin_cal_frd);
    Zin_il_pll_avg_cal_frd = (Ypll_il_frd+Yin_cal_frd)\T_frd;
    Zin_il_pll_avg_cal = (Ypll_il+Yin_cal)\T;
    Ti_frd = Gid_cal_frd*Gdel_frd*(-Gdei_frd + Gci_frd)*Ki_frd;

%     figure(4)
%     bode(Zin_il_pll_avg_cal_frd,Bode_O)
% %     bode(Zin_il_pll_avg_cal_frd,Zo_il_pll_avg_sim,Bode_O)
%     hold on
%     Bode_Darklines(3)
%     Zout_VSI_il_pll_qq=Zin_il_pll_avg_cal_frd(2,2);
%     figure(41)
%     bode(Zout_VSI_il_pll_qq,Bode_O)
% %     bode(Zin_il_pll_avg_cal_frd,Zo_il_pll_avg_sim,Bode_O)
%     hold on
%     Bode_Darklines(3)
%     figure(42)
%     pzmap(Zo_il_pll_avg_sim(2,2))
%     hold on
%     
%     figure(5)
%     bode(Zin_il_pll_avg_cal_frd,Z_RLC_dq,Zo_il_pll_avg_sim,Bode_O)
% %     bode(Zin_il_pll_avg_cal_frd,Zo_il_pll_avg_sim,Bode_O)
%     hold on
%     Bode_Darklines(3)
    figure(6)
    Zout_VSI_il_pll_qq = Zin_il_pll_avg_cal_frd(2,2);
    
    Z_RLC_qq = Z_RLC_dq(2,2);
    bode(Zout_VSI_il_pll_qq,Z_RLC_qq,Bode_O)
    hold on
    Bode_Darklines(3)
    figure(611)
%     Zout_VSI_qq =Vsd_s/Id+s^2/(Id*Kp1*s+Id*Ki1);
    Zout_VSI_qq = (s/tf_pll+Vsd_s)/Id;
    Zout_VSI_qq = (s+Vsd_s*tf_pll)/(Id*tf_pll);
    N=0.015;
    Zout_VSI_qq_N = (s/tf_pll+Vsd_s)/(Id*(N*s+1));
%     Zout_VSI_qq_N = 1/(Id*((1-N*tf_pll*Vsd_s)*Gpll+N*tf_pll));
    Zout_VSI_qq_v1 = (s/tf_pll+Vsd_s)/(Dd/Gci(2,2)/Vdc+Id);
    bode(Zout_VSI_il_pll_qq,Zout_VSI_qq,Zout_VSI_qq_N,Bode_O)
    Bode_Darklines(3)
    hold on
    
    figure(612)
    bode(Zout_VSI_il_pll_qq,Zout_VSI_qq_v1,Bode_O)
    Bode_Darklines(3)
    hold on
    
%     figure(61)
%     nyquist(Z_RLC_qq/Zout_VSI_il_pll_qq);
%     hold on
%     figure(62)
%     nyquist(Zout_VSI_il_pll_qq/Z_RLC_qq);
%     hold on
%     figure(63)
%     pzmap(Zo_il_pll_avg_sim(2,2)/Z_RLC_qq)
%     hold on
% %     figure(62)
% %     pzmap(Z_RLC_qq/Zout_VSI_il_pll_qq)
% %     hold on
%     
% %     Zo_il_pll_avg_sim=[Zo_il_pll_avg_sim(1,1) 0; 0 Zo_il_pll_avg_sim(2,2)];
% %     Z_RLC_dq=[Z_RLC_dq(1,1) 0; 0 Z_RLC_dq(2,2)];
%     LL=(Zo_il_pll_avg_sim/Z_RLC_dq);
%     n=1e3;
%     f=logspace(-4,4,n);
%     w=f*2*pi;
%     % For negative frequencies turn on
% %     if(1)
% %         w=[-fliplr(w) w];
% %     %     n=n*2;
% %     end
%     % L(s) frequency response
%     Lresp=freqresp(LL,w);
%     for k=1:length(w)
%         Leigenvalues(:,k)=eig(Lresp(:,:,k));
%     end
% 
%     figure(7);
%     plot(1./Leigenvalues(1,:),linestyle(i),'LineWidth',linewidth)
%     hold on
%     plot(1./Leigenvalues(2,:),strcat(':',linestyle(i)),'LineWidth',linewidth)
%     grid on
%     axisloci=axis;
% %     legend({'{\it\lambda}_{1}','{\it\lambda}_{2}','Ldd'},'Fontsize',fontsize,'FontWeight','bold')
%     plot(-1,0,'r+','LineWidth',linewidth)
% %     hold off
%     title('Characteristic Loci of L','Fontsize',fontsize,'FontWeight','bold')
%     grid on
%     set(gca,'FontSize',fontsize);
%     
%     h=ezplot('x^2+y^2=1');
%     set(h,'color','r');
%     
%     figure(8)
%     pzmap(LL)
%     hold on
end
