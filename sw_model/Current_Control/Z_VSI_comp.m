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

load Z_VSI_kp_1_ki_2_N_0.mat
Z_VSI_N_0=Z.ZDQL;
load Z_VSI_kp_1_ki_2_N_005_v1.mat
Z_VSI_N_005=Z.ZDQL;
load Z_VSI_kp_1_ki_2_N_015_v1.mat
Z_VSI_N_015=Z.ZDQL;
figure(1)
bode(Z_VSI_N_0/3,Z_VSI_N_005/3,Z_VSI_N_015/3,Bode_O)

load Z_VSI_kp_1_7_ki_158_N_0.mat
Z_VSI_N_0=Z.ZDQL;
load Z_VSI_kp_1_7_ki_158_N_005.mat
Z_VSI_N_005=Z.ZDQL;
load Z_VSI_kp_1_7_ki_158_N_0065.mat
Z_VSI_N_0065=Z.ZDQL;
figure(2)
bode(Z_VSI_N_0/3,Z_VSI_N_005/3,Z_VSI_N_0065/3,Bode_O)