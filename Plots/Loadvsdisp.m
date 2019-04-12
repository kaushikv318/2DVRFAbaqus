clc
clear all
%close all

%%
fid0 = load('../0_interface/RF_0int.rpt');
fid1 = load('../1_interface/RF_1int.rpt');
fid2 = load('../2_interface/RF_2int.rpt');
fid3 = load('../3_interface/RF_3int.rpt');
fid4 = load('../4_interface/RF_4int.rpt');

rf_0 = -fid0(:,2);
u_0 = -fid0(:,3);
rf_1 = -fid1(:,2);
u_1 = -fid1(:,3);
rf_2 = -fid2(:,2);
u_2 = -fid2(:,3);
rf_3 = -fid3(:,2);
u_3 = -fid3(:,3);
rf_4 = -fid4(:,2);
u_4 = -fid4(:,3);

fps =24;


f1=figure(1)
plot(u_0,rf_0,'k',u_1,rf_1,'r',u_2,rf_2,'b',u_3,rf_3,'m',u_4,rf_4,'c','linewidth',1,'markersize',2)
xlabel('$\Delta$(mm)','Fontsize',14,'Interpreter','latex')
ylabel('$P$(N)','Fontsize',14,'Interpreter','latex')
title('Load vs Displacement','Interpreter','latex','Fontsize',14) 
set(gca,'Fontsize',14,'Position',[0.13 0.15 0.775 0.76]);
h1 = legend('0 interface','1 interface','2 interface','3 interface','4 interface');
set(h1,'Box','Off','Fontsize',8,'FontName','Helvetica','Location','Best');
axis([0 0.2 0 40])
pbaspect([1.5 1.2 1])
publish_fig(1,fullfile('RF'),f1,'','','',2.6,8)

hold on

