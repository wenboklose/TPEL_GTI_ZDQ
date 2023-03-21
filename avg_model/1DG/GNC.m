close all
clear all
%% Bode plot options
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e3]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth=3; fontsize=18;
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='m';
linestyle5='k'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3];
s=tf([1 0],[0 1]);
%% passive load
RLoad = 10*tf([0 1],[0 1]);
CLoad = 250e-6;
Z_RLoad_dq = [RLoad 0;0 RLoad];
Y_CLoad_dq = [CLoad*s -60*2*pi*CLoad;60*2*pi*CLoad CLoad*s];
Z_RC_dq = 1/(1/Z_RLoad_dq+Y_CLoad_dq);

%% RL grid
R_grid = .2;
L_grid = 2e-3;
Z_grid = R_grid + s*L_grid;
Z_grid_dq = [Z_grid -1*60*2*pi*L_grid; 1*60*2*pi*L_grid Z_grid];

load Zout_VSI_190_s.mat %Z_out_VSI_u.mat
Z_VSI=Zin_il_pll_avg_sim;

Z_S = 1/(1/Z_VSI);
Z_L = 1/(1/Z_RC_dq+1/Z_grid_dq);

% %% R grid
% R_grid = 1.2;
% Z_grid = R_grid;
% Z_grid_dq = JF_DQFromABC(Z_grid,60*2*pi);
% 
% load Zout_VSI_170_u_R_grid.mat %Z_out_VSI_u.mat
% Z_VSI=Zin_il_pll_avg_sim;
% load Zin_AFE_170_R_grid.mat %Z_in_AFE_u_v300vdc.mat %Z_in_AFE_s_v1.mat %
% Z_AFE=Zin_vl_pll_avg_sim;
% 
% Z_S = Z_VSI;
% 
% Z_L = 1/(1/Z_AFE+1/Z_RC_dq+1/Z_grid_dq);

%     LL=(Z_L)/Z_S;
    LL=Z_S/Z_L;

% plot GNC for using RR and L
    n=4e3;
    f=logspace(-8,8,n);
%     f=1:1:1e8;
    w=f*2*pi;
    % For negative frequencies turn on
%     if(1)
%         w=[-fliplr(w) w];
%     %     n=n*2;
%     end
    % L(s) frequency response
    Lresp=freqresp(LL,w);
    for k=1:length(w)
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
        Leigenvalues(:,k)=eig(Lresp(:,:,k));
%         leigenvalues(:,k)=eig(rr_res(:,:,k));
    end
%     [Leigenvalues]=sortloci(Leigenvalues);
    figure(7);
    plot(Leigenvalues(1,:),linestyle(1),'LineWidth',linewidth)
    hold on
    plot(Leigenvalues(2,:),strcat(':',linestyle(2)),'LineWidth',linewidth)
    
    figure(7);
    plot(1./Leigenvalues(1,:),linestyle(1),'LineWidth',linewidth)
    hold on
    plot(1./Leigenvalues(2,:),strcat(':',linestyle(2)),'LineWidth',linewidth)
    hold on
    grid on
    axisloci=axis;
    legend({'{\it\lambda}_{1}','{\it\lambda}_{2}','Ldd'},'Fontsize',fontsize,'FontWeight','bold')
    plot(-1,0,'r+','LineWidth',linewidth)
%     hold off
    title('Characteristic Loci of L','Fontsize',fontsize,'FontWeight','bold')
    grid on
    set(gca,'FontSize',fontsize);
    
    h=ezplot('x^2+y^2=1');
    set(h,'color','r');
    
    figure(8)
    pzmap(LL)
    hold on
    
    figure(9)
    bode(Z_S,Z_L,Bode_O)
    Bode_Darklines(3)
    hold on
    Z_S_qq=Z_S(2,2);
    Z_L_qq=Z_L(2,2);
    figure(10)
    bode(Z_S_qq,Z_L_qq,Bode_O)
    Bode_Darklines(3)
    hold on
    figure(11)
    nyquist(Z_L_qq/Z_S_qq,{1e-2,1e4})
    Bode_Darklines(3)
    hold on
    figure(12)
    nyquist(Z_S_qq/Z_L_qq,{1e-2,1e8})
    Bode_Darklines(3)
    hold on
    figure(13)
    pzmap(Z_S_qq/Z_L_qq)
    hold on
%     
%     figure(1)
%     bode(Z_VSI,Bode_O)
%     Bode_Darklines(3)
% 
%     figure(2)
%     bode(Z_AFE,Bode_O)
%     Bode_Darklines(3)
%     hold on
%     figure(3)
%     bode(Z_RC_dq,Bode_O)
%     Bode_Darklines(3)
 %% calculate the ACindex   
%     Z_VSI_s = [Z_VSI(1,1) tf([0 0],[0 1]);tf([0 0],[0 1]) Z_VSI(2,2)];
%     Z_L_s = [Z_L(1,1) tf([0 0],[0 1]);tf([0 0],[0 1]) Z_L(2,2)];
%     ll = Z_L_s/Z_VSI;
%     n=1e2;
%     f=logspace(-4,4,n);
%     % f=40:.05:1e3;
%     w=f*2*pi;
%     % L(s) frequency response
%     llresp=freqresp(ll,w);
%     RRdd=frd(squeeze(llresp(1,1,1:n)),w);
%     RRdq=frd(squeeze(llresp(1,2,1:n)),w);
%     RRqd=frd(squeeze(llresp(2,1,1:n)),w);
%     RRqq=frd(squeeze(llresp(2,2,1:n)),w);
%     % calculate ACindex defined by Rolando
%     ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);
%     figure(4)
%     bode(ACindex,Bode_O)
%     Bode_Darklines(3)
%     
%     
%     LLresp=freqresp(1/LL,w);
%     RRdd=frd(squeeze(llresp(1,1,1:n)),w);
%     RRdq=frd(squeeze(llresp(1,2,1:n)),w);
%     RRqd=frd(squeeze(llresp(2,1,1:n)),w);
%     RRqq=frd(squeeze(llresp(2,2,1:n)),w);
%     LL_v1=[RRdd RRdq;RRqd RRqq];
%     figure(5)
%     bode(LL_v1*LL_v1',LL_v1'*LL_v1,Bode_O)
%     Bode_Darklines(3)