%% plot time domain simulation results
close all
time = f_vsi(1:123247,1);
fvsi = f_vsi(1:123247,2);
% fafe = fvsi_afe(1:223247,3);
va   = Vpcc(1:123247,2);
vb   = Vpcc(1:123247,3);
vc   = Vpcc(1:123247,4);

ia   = IVSI(1:123247,2);
ib   = IVSI(1:123247,3);
ic   = IVSI(1:123247,4);

figHandle=figure(1);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(time,fvsi,'r')
% hold on
% plot(time,fafe,'k')
grid on
xlabel('Time (s)','FontSize',20)
ylabel('Frequency (Hz)','FontSize',20)
xlim([9.9 10.1])
ylim([40 80])
% title('Ac input current of AFE','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')

figHandle=figure(2);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(time,va,'r')
hold on
plot(time,vb,'k')
plot(time,vc,'b')
grid on
xlabel('Time (s)','FontSize',20)
ylabel('Voltage (V)','FontSize',20)
xlim([9.9 10.1])
ylim([-450 450])
% title('Ac input current of AFE','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')

figHandle=figure(3);
set(figHandle, 'Position', [50, 10, 1000, 800]);
plot(time,ia,'r')
hold on
plot(time,ib,'k')
plot(time,ic,'b')
grid on
xlabel('Time (s)','FontSize',20)
ylabel('Current (A)','FontSize',20)
xlim([9.9 10.1])
ylim([-450 450])
% title('Ac input current of AFE','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')

% figure(2)
% plot(vll_vdc.time,vll_vdc.signals(2).values)
% grid on
% xlim([0 1])
% ylim([0 3])
% xlabel('Time (s)','FontSize',20)
% ylabel('Voltage (V)','FontSize',20)
% title('Dc interface voltage','FontSize',20)
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3);
% set(gca,'FontSize',20)
% set(gca,'Fontname','Times New Roman')