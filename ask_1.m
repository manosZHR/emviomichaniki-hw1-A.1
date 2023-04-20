clear; close all; clc

Cinf=0.04; 
C_critical=0.01;
Deff=0.01; 
r_origin=0; r_last=400;
to=3600;
k=to*Deff/(r_last^2);
par = [Cinf k];

nn_rn=r_last;
rn_origin=0; rn_last=1; drn = (rn_last-rn_origin)/nn_rn;
rn=linspace(rn_origin,rn_last,nn_rn);

%% TASK A
c0=zeros(1,nn_rn);
c0(nn_rn)=Cinf;

days=7;
tau_origin=0; tau_last=days*24;
steps=tau_last;
tau_span=linspace(tau_origin,tau_last,steps);


[tau_i,ci]=ode89(@(tau,c)odeA1(tau,c,par,nn_rn),tau_span,c0);

figure(1)
plot(r_last.*rn,ci)
hold on
yl=yline(C_critical,'--k','critical concentration','LineWidth',2);
yl.LabelHorizontalAlignment = 'left';
hold off
grid on; box on; axis tight
xlabel('r (μm)'); ylabel('concenctration (μg/mL)');

figure(2)
for i=10:10:tau_last
    plot(r_last.*rn,ci(i,:),'DisplayName',strcat('t=',num2str(i)))
    hold on
end
yl=yline(C_critical,'--k','critical concentration','LineWidth',2);
yl.LabelHorizontalAlignment = 'center';
hold off
grid on; box on; axis tight
xlabel('r (μm)'); ylabel('concenctration (μg/mL)'); legend('show','Location','west')

figure(3)
M = readmatrix('comsol_6_24_168.txt');
plot(M(1:10:301,1),M(1:10:301,2),'or','MarkerSize',5)
hold on
plot(M(302:10:602,1),M(302:10:602,2),'og','MarkerSize',5)
plot(M(603:10:903,1),M(603:10:903,2),'ob','MarkerSize',5)

plot(r_last.*rn,ci(6,:),'r')
plot(r_last.*rn,ci(24,:),'g')
plot(r_last.*rn,ci(168,:),'b')

yl=yline(C_critical,'--k','critical concentration','LineWidth',2);
yl.LabelHorizontalAlignment = 'left';
hold off
grid on; box on; axis tight
xlabel('r (μm)'); ylabel('concenctration (μg/mL)');

lgd=cell(6,1);
lgd{1}='6 h (comsol)';
lgd{2}='24 h (comsol)';
lgd{3}='1 w (comsol)';
lgd{4}='6 h (matlab)';
lgd{5}='24 h (matlab)';
lgd{6}='1 w (matlab)';
legend(lgd,'Location','northwest')

%% TASK B

c0=ci(end,:);
c0(nn_rn)=0;

[tau_i,ci]=ode89(@(tau,c)odeA1(tau,c,par,nn_rn),tau_span,c0);

figure(4)
plot(r_last.*rn,ci)
hold on
yl=yline(C_critical,'--k','critical concentration','LineWidth',2);
yl.LabelHorizontalAlignment = 'left';
hold off
grid on; box on; axis tight
xlabel('r (μm)'); ylabel('concenctration (μg/mL)');

figure(5)
for i=10:10:tau_last
    plot(r_last.*rn,ci(i,:),'DisplayName',strcat('t=',num2str(i)))
    hold on
end
yl=yline(C_critical,'--k','critical concentration','LineWidth',2);
yl.LabelHorizontalAlignment = 'center';
hold off
grid on; box on; axis tight
xlabel('r (μm)'); ylabel('concenctration (μg/mL)'); legend('show','Location','west')

figure(6)
plot(r_last.*rn,ci(6,:),'r')
hold on
plot(r_last.*rn,ci(24,:),'g')
plot(r_last.*rn,ci(168,:),'b')

yl=yline(C_critical,'--k','critical concentration','LineWidth',2);
yl.LabelHorizontalAlignment = 'left';
hold off
grid on; box on; axis tight
xlabel('r (μm)'); ylabel('concenctration (μg/mL)');

lgd=cell(3,1);
lgd{1}='6 h';
lgd{2}='24 h';
lgd{3}='1 w';
legend(lgd,'Location','northwest')
