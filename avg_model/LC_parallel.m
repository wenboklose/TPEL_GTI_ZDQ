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
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='m';
linestyle5='k'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3];
s=tf([1 0],[0 1]);
L=1000e-6;C=50e-6;
Z_LC = 1/(1/(L*s)+1/(1/(C*s)));
Z_LC_dq = JF_DQFromABC(Z_LC,400*2*pi);

figure(1)
bode(Z_LC_dq,Bode_O)