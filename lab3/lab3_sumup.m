clc
close all
clear all
%*****************************prelab****************************%
%The component in the transfer function [-a1, -a2, b1, b2]
% system using dataset SysIdenData.mat
% est_val = [1.6094, -0.6176, 0.001026, 0.01291];

% system using dataset SysIdenData_StudentVersion.mat
% est_val = [1.3551, -0.3707, 0.0095 0.0364];

% system using dataset SysIdenData_lab4.mat
est_val = [1.2512, -0.2724, 0.0147, 0.0412];

% recreate of transfer function and state space model
Ts = 0.75;
G = [0, 1; est_val(2), est_val(1)];
H = [0;1];
C = [est_val(4), est_val(3)];
D = 0;

sys1 = tf([est_val(3), est_val(4)], [1, -est_val(1), -est_val(2)], Ts, 'Variable', 'z^-1')
sys2 = ss(G, H, C, D, Ts)

% find the eigenvalue for the polynomial in the denominator term
% pole(sys1)
% roots([1, -est_val(1), -est_val(2)])
eig(G)

% find the zero for the polynomial in the numerator term
[z, g] = zero(sys1)
% non minimum phase system( the zero is outside the unit circle)

% determine the rank of controllability matrix [H GH] full rank this case->controllable
str = ['The rank of [H,G*H] is ',num2str(rank([H, G*H]))];
disp(str)

% reachable as well for continuous and dicrete model
si = size([H, G*H]);
if (rank([H, G*H])== min(si(1), si(2)))
    disp('The system is reachable as well as controllable');
end

% determine the rank of observability matrix [C; CG]

str = ['The rank of [C,C*G] is ',num2str(rank([C, C*G]))];
disp(str)
si = size([C, C*G]);
if (rank([C, C*G])== min(si(1), si(2)))
    disp('The system is observable as well as detectable');
end

G_D = G';
H_D = C';
C_D = H';
D_D = D';
rank([H_D, G_D*H_D])
rank([C_D; C_D*G_D])

%*****************************part 1****************************%

% part a
% set x1 = 10
x = [0 0.3]'
% the abs(u(k)max) is 0.575
% for deadbeat control
p0 = 0;
p1 = 0;
L_db = acker(G_D, H_D, [p0, p1])
W_C = [H_D, G_D*H_D]
L_db2 = [0, 1]*inv(W_C)*G_D*G_D
%for non-deatbeat control

% p0 = 0.6315;
% p1 = 0.9775;
p0 = 0.4;
p1 = 0.9;
L_ndb = acker(G_D, H_D, [p0, p1])
sim('lab3_e1');

%*****************************part 2****************************%
% calculation of DC gain

p0 = 0.9;
p1 = 0.92;
L_ndb = acker(G_D, H_D, [p0, p1])
sys3 = ss(G_D-H_D*L_ndb, H_D, C_D-D_D*L_ndb, D_D, Ts)
Gcl = dcgain(sys3)
sim('lab3_e2');

figure(1)
subplot(2,1,1)
plot(ScopeData(:,1),ScopeData(:,4),'r');hold on;
stairs(ScopeData(:,1),ScopeData(:,2),'g');hold off;
xlim([0,600]);ylim([-0.8,0.8]);grid on;
title({'Set-point control results:simulation','Output signal'});
xlabel({'Time(sec)','(a)'});ylabel({'Offset-Free';'Water Level(V)'});
legend('Simulated Output','Reference Output','location','best');

subplot(2,1,2)
plot(ScopeData(:,1),ScopeData(:,3),'b');
grid on;
title('Offset-Free Control Input signal');
xlabel({'Time(sec),(b)'});
ylabel({'Ofefset-Free';'Pump Volatge(V)'});
legend({'Simulated control Output'},'location','best');
xlim([0,600]);ylim([-0.5,0.5]);

%*****************************part 3****************************%
% p0 = 0.9;
% p1 = 0.9;
% L_ndb = acker(G_D, H_D, [p0, p1])
p0 = 0;
p1 = 0;
K = acker(G, H, [p0 p1])
K = K';
sim('lab3_e3');
figure(2)
plot( ScopeData1.time,ScopeData1.signals.values(:,1),'r');hold on
stairs( ScopeData1.time,ScopeData1.signals.values(:,2),'b');hold off
grid on;
title('\bf State Estimation Error');
xlabel('Time(sec),(C)');
ylabel({'\bf Estimation Error'});
legend({'x_{1}(k)-x_{1}^*(k)','x_{2}(k)-x_{2}^*(k)'});
xlim([0,3]);ylim([0,0.4]);

%*****************************post lab****************************%
%*****************************part a******************************%
W = [0.5, 0]';
% dcgain changes due to the introduction of disturbance
% sys4 = ss(G_D-H_D*L_ndb, W-H_D, C_D, D_D, Ts);
sys4 = ss(G_D-H_D*L_ndb, H_D, C_D, D_D, Ts);
Gcl = dcgain(sys4)
sim('lab3_postlab');
figure(3)
plot(ScopeData2(:,1),ScopeData2(:,4),'r');hold on;
stairs(ScopeData2(:,1),ScopeData2(:,3),'g');hold off;
xlim([0,600]);grid on;
title({'Set-point control results:simulation','Output signal'});
xlabel({'Time(sec)','(a)'});ylabel({'Offset-Free';'Water Level(V)'});
legend('Simulated Output','Reference Output','location','best');

%*****************************part c******************************%
% TF=(C(Iz-(G-H*L))^(-1)*(W-H*Lv) TF=0 | z = 1
temp = G_D-H_D*L_ndb;
L_v=(1-temp(2,2))/(H_D(1)*(1-temp(2,2))+H_D(2)*temp(1,2));
% L_v = 0.2771;
sys5 = ss(G_D-H_D*L_ndb, W-H_D*L_v, C_D, D_D, Ts);
Gcl = dcgain(sys5);
G_a = [G_D W; 0 0 1];
H_a = [H_D; 0];
C_a = [C_D 0];

% K_D = G_a^3*inv([C_a; C_a*G_a; C_a*G_a*G_a])*[0; 0; 1]
% sim('lab3_postlab_bonus');


est_val_err = est_val+(rand(1,4)*0.05-0.1);
G = [0, 1; est_val_err(2), est_val_err(1)];
H = [0;1];
C = [est_val_err(4), est_val_err(3)];
D = 0;
G_D = G';
H_D = C';
C_D = H';
D_D = D';
% sys3 = ss(G_D-H_D*L_ndb, H_D, C_D-D_D*L_ndb, D_D, Ts)
% Gcl = dcgain(sys3)
sim('lab3_postlab2');