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
R=1/(1/10+1/10);C=50e-6;
Z_RC = R/(R*C*s+1);
Z_RC_dq = JF_DQFromABC(Z_RC,60*2*pi);

load Z_VSI_u.mat
Z_VSI=Zin_il_pll_avg_sim;
load Z_AFE_u.mat
Z_AFE=Zin_il_pll_avg_sim;
Z_L = 1/(1/Z_AFE+1/Z_RC_dq);
% Z_VSI = [Z_VSI(1,1) 0; 0 Z_VSI(2,2)];
% Z_L = [Z_L(1,1) 0; 0 Z_L(2,2)];
    LL=Z_VSI/(Z_L);
    LL=(Z_L)/Z_VSI;
%     L=Z_VSI/Z_AFE;
% plot GNC for using RR and L
    n=2e3;
    f=logspace(-4,4,n);
    % f=40:.05:1e3;
    w=f*2*pi;
    % For negative frequencies turn on
    if(1)
        w=[-fliplr(w) w];
    %     n=n*2;
    end
    % L(s) frequency response
    Lresp=freqresp(LL,w);
%     lresp=freqresp(l,w);
%     Lddresp=freqresp(Ldd,w);
%     Lqqresp=freqresp(Lqq,w);
    % Eigenvalues and Sorting
%     RR_res=[RR(1,1).ResponseData RR(1,2).ResponseData;RR(2,1).ResponseData RR(2,2).ResponseData];
%     rr_res=[rr(1,1).ResponseData rr(1,2).ResponseData;rr(2,1).ResponseData rr(2,2).ResponseData];
    for k=1:length(w)
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
        Leigenvalues(:,k)=eig(Lresp(:,:,k));
%         leigenvalues(:,k)=eig(rr_res(:,:,k));
    end
    % [Leigenvalues]=sortloci(Leigenvalues);
    % Characteristic Loci
%     figure(5)
%     %clf
%     hold on
%     plot(Lddeigenvalues(1,:),linestyle1,'LineWidth',linewidth)
%     % hold on
%     % plot(Lddeigenvalues(2,:),linestyle2,'LineWidth',linewidth)
%     legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
%     plot(-1,0,'r+','LineWidth',linewidth)
%     hold off
%     title('Characteristic Loci of Ldd','Fontsize',fontsize,'FontWeight','bold')
%     grid
%     set(gca,'FontSize',fontsize);
% 
%     figure(6)
%     %clf
%     hold on
%     plot(Leigenvalues(1,:),'LineWidth',linewidth)
%     hold on
%     % plot(Leigenvalues(2,:),linestyle2,'LineWidth',linewidth)
%     legend({'{\it\lambda}_{1}','{\it\lambda}_{2}'},'Fontsize',fontsize,'FontWeight','bold')
%     plot(-1,0,'r+','LineWidth',linewidth)
% %     hold off
%     title('Characteristic Loci of L','Fontsize',fontsize,'FontWeight','bold')
%     grid
%     set(gca,'FontSize',fontsize);


    figure(7);
    plot(Leigenvalues(1,:),linestyle(1),'LineWidth',linewidth)
    hold on
    plot(Leigenvalues(2,:),strcat(':',linestyle(2)),'LineWidth',linewidth)
    % plot(Lddeigenvalues(1,:),linestyle3,'LineWidth',linewidth)
    % plot(Lqqeigenvalues(1,:),linestyle4,'LineWidth',linewidth)
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
    bode(Z_VSI,Z_L,Bode_O)