%% plot time domain simulation results
% tf_pll.time = 0:0.1:1;
% tf_pll.signals(1).values = 0:0.1:1;
% tf_pll.signals(2).values = 0:0.1:1;
figure(1)
plot(f_pll.time,f_pll.signals(1).values)
grid on
xlabel('Time (s)','FontSize',20)
ylabel('Current (V)','FontSize',20)
xlim([0 0.6])
ylim([250 500])
title('VSI PLL Frequency','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')
figure(2)
plot(va_ia.time,va_ia.signals(1).values)
grid on
xlim([0 0.6])
ylim([-100 100])
xlabel('Time (s)','FontSize',20)
ylabel('Voltage (V)','FontSize',20)
title('Ac Interface Voltage','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')
figure(3)
plot(va_ia.time,va_ia.signals(2).values)
grid on
xlim([0 0.6])
ylim([-12 12])
xlabel('Time (s)','FontSize',20)
ylabel('Current (A)','FontSize',20)
title('VSI Output Current','FontSize',20)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3);
set(gca,'FontSize',20)
set(gca,'Fontname','Times New Roman')