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
    Bode_O.XLim={[1 1e4]};
    Bode_O.XLimMode={'manual'};
    Bode_O.FreqUnits='Hz';
    Bode_O.PhaseWrapping='on';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    linewidth=3; fontsize=18;
    linestyle1='k'; linestyle2='r';
    linestyle3='r'; linestyle4='m';
    linestyle5='k'; linestyle6='c';

    linestyle=[linestyle1; linestyle2; linestyle3];
    s=tf([1 0],[0 1]);
%% RL grid
    R_grid = .2;
    L_grid = 2e-3;
    Z_grid = R_grid + s*L_grid;
    Z_grid_dq = [Z_grid -1*60*2*pi*L_grid; 1*60*2*pi*L_grid Z_grid];

%% passive load
    RLoad = 10*tf([0 1],[0 1]);
    CLoad = 250e-6;
    Z_RLoad_dq = [RLoad 0;0 RLoad];
    Y_CLoad_dq = [CLoad*s -60*2*pi*CLoad;60*2*pi*CLoad CLoad*s];
    Z_RC_dq = 1/(1/Z_RLoad_dq+Y_CLoad_dq);
    load Zout_VSI_190_s.mat %Z_out_VSI_u.mat
    Z_VSI=Zin_il_pll_avg_sim;

    Z_L = 1/(1/Z_RC_dq+1/Z_grid_dq);
    Z_Load_dd = Z_L(1,1);
    Z_Load_dq = Z_L(1,2);
    Z_Load_qd = Z_L(2,1);
    Z_Load_qq = Z_L(2,2);
    Z_Source_s = Z_VSI;
    Z_Source_dd_s = Z_Source_s(1,1);
    Z_Source_dq_s = Z_Source_s(1,2);
    Z_Source_qd_s = Z_Source_s(2,1);
    Z_Source_qq_s = Z_Source_s(2,2);
%     Z_Source_s2 = Z_VSI_s2;
%     Z_Source_dd_s2 = Z_Source_s2(1,1);
%     Z_Source_dq_s2 = Z_Source_s2(1,2);
%     Z_Source_qd_s2 = Z_Source_s2(2,1);
%     Z_Source_qq_s2 = Z_Source_s2(2,2);
    
    load Zout_VSI_190_u.mat %Z_out_VSI_u.mat
    Z_VSI=Zin_il_pll_avg_sim;
    
    Z_Source_u = Z_VSI;
    Z_Source_dd_u = Z_Source_u(1,1);
    Z_Source_dq_u = Z_Source_u(1,2);
    Z_Source_qd_u = Z_Source_u(2,1);
    Z_Source_qq_u = Z_Source_u(2,2);
    
    figHandle=figure(1);
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Z_Source_dd_s,Z_Load_dd,Bode_O)
    hold on
    bode(Z_Source_dd_u,Bode_O)
    Bode_Darklines(3)
        
    figHandle=figure(2);
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Z_Source_dq_s,Z_Load_dq,Bode_O)
    hold on
    bode(Z_Source_dq_u,Bode_O)
    Bode_Darklines(3)
    
    figHandle=figure(3);
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Z_Source_qd_s,Z_Load_qd,Bode_O)
    hold on
    bode(Z_Source_qd_u,Bode_O)
    Bode_Darklines(3)

    figHandle=figure(4);
    set(figHandle, 'Position', [50, 10, 1000, 800]);
    bode(Z_Source_qq_s,Z_Load_qq,Bode_O)
    hold on
%     bode(Z_Source_qq_s2,Z_Load_qq_s,Bode_O)
    bode(Z_Source_qq_u,Bode_O)
    Bode_Darklines(3)