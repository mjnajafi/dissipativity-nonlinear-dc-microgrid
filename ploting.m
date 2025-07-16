close all;
fontsize=35;
axiswidth=2;
graphwidth=2.5;
min_t=0;
max_t=10;

%% Voltage and per-unit current at the same time
figure(1)
subplot(211)

V1=plot(out.V(:,1), out.V(:,2), '-r'); 
hold on;
V2=plot(out.V(:,1), out.V(:,3), '--g'); 
V3=plot(out.V(:,1), out.V(:,4), ':b'); 
V4=plot(out.V(:,1), out.V(:,5), '-.m'); 
V5=plot(out.V(:,1), out.V(:,6), '--c'); 
% V6=plot(out.V(:,1), out.V(:,7), ':y'); 
% V7=plot(out.V(:,1), out.V(:,8), '-k'); 
y=ylabel('\bf Voltages(V)');
% x=xlabel('\bf Time(sec)'); 
% l=legend({'\bf DG1','\bf DG2','\bf DG3','\bf DG4','\bf DG5','\bf DG6'},...
%     'Location','southeast','Orientation','horizontal');
grid on
title('\bf (a)')
set(V1,'LineWidth',graphwidth);
set(V2,'LineWidth',graphwidth);
set(V3,'LineWidth',graphwidth);
set(V4,'LineWidth',graphwidth);
set(V5,'LineWidth',graphwidth);
% set(V6,'LineWidth',graphwidth);
% set(V7,'LineWidth',graphwidth);
% set(x,'FontSize',fontsize,'FontName','Times New Roman');
set(y,'FontSize',fontsize,'FontName','Times New Roman');
% set(l,'FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',fontsize, 'LineWidth', axiswidth,'fontweight','bold',...
    'xlim',[min_t max_t],'ylim',[42 52],'FontName','Times New Roman','xticklabel',[]);

subplot(212) 

p1=plot(out.P(:,1), out.P(:,2), '-r'); 
hold on;
p2=plot(out.P(:,1), out.P(:,3), '--g'); 
p3=plot(out.P(:,1), out.P(:,4), ':b'); 
p4=plot(out.P(:,1), out.P(:,5), '-.m'); 
% p5=plot(out.P(:,1), out.P(:,6), '--c'); 
% p6=plot(out.P(:,1), out.P(:,7), ':y'); 
y=ylabel('\bf Currents(p.u.)');
x=xlabel('\bf Time(sec)'); 
l=legend({'\bf DG1','\bf DG2','\bf DG3','\bf DG4','\bf DG5','\bf DG6'},...
    'Location','southeast','Orientation','horizontal');
grid on
title('\bf (b)')
set(p1,'LineWidth',graphwidth);
set(p2,'LineWidth',graphwidth);
set(p3,'LineWidth',graphwidth);
set(p4,'LineWidth',graphwidth);
% set(p5,'LineWidth',graphwidth);
% set(p6,'LineWidth',graphwidth);
set(x,'FontSize',fontsize,'FontName','Times New Roman');
set(y,'FontSize',fontsize,'FontName','Times New Roman');
set(l,'FontSize',30,'FontName','Times New Roman');
set(gca,'FontSize',fontsize, 'LineWidth', axiswidth,'fontweight','bold',...
    'xlim',[min_t max_t],'ylim',[0 1],'FontName','Times New Roman');

%% Average Voltage
% figure(1)
% AveVoltage=plot(out.Ave(:,1),out.Ave(:,2),'-');
% hold on
% RefVoltage=plot(out.V(:,1),out.V(:,6));
% y=ylabel('\bf Average Voltage(V)');
% x=xlabel('\bf Time(sec)'); 
% grid on
% set(AveVoltage,'LineWidth',graphwidth);
% set(RefVoltage,'LineWidth',1);
% set(x,'FontSize',fontsize,'FontName','Times New Roman');
% set(y,'FontSize',fontsize,'FontName','Times New Roman');
% set(gca,'FontSize',fontsize, 'LineWidth', axiswidth,'fontweight','bold',...
%     'xlim',[min_t max_t],'ylim',[0 60],'FontName','Times New Roman');
% hold on


%% Average Voltage
% figure(5)
% Average=plot(out.Avee(:,1),out.Avee(:,2),'-');
% y=ylabel('\bf Average Voltage(V)');
% x=xlabel('\bf Time(sec)'); 
% grid on
% set(Average,'LineWidth',graphwidth);
% set(x,'FontSize',fontsize,'FontName','Times New Roman');
% set(y,'FontSize',fontsize,'FontName','Times New Roman');
% set(x,'FontSize',fontsize,'FontName','Times New Roman');
% set(gca,'FontSize',fontsize, 'LineWidth', axiswidth,'fontweight','bold',...
%     'xlim',[min_t max_t],'ylim',[3 66],'FontName','Times New Roman');
% hold on
