% system identification
clc;
close all;
clear all;
load('SysIdenData_lab5.mat');
y_actm = LogData.signals(1).values(:,1);
y_act = LogData.signals(1).values(:,2);
u_act = LogData.signals(2).values;
t = LogData.time;
Ts = t(2)-t(1);

% make trucation here
trunc_time = 139.5;
cycle_time = 600;
y_actm = y_actm(trunc_time/Ts:2*cycle_time/Ts,:);
y_act = y_act(trunc_time/Ts:2*cycle_time/Ts,:);
u_act = u_act(trunc_time/Ts:2*cycle_time/Ts,:);
t = t(trunc_time/Ts:2*cycle_time/Ts,:);

i = 1;
ave = 0;
while(u_act(i+1) == u_act(i))
%     ave = mean(y_act(1:i,:));
    i = i+1;
end
% y_offset
ave=mean(y_act(1:i))
y = y_act-ave;
u = u_act-u_act(1);

start = 3;
N = ceil(size(y, 1)/2);

Matrix(:, 1) = y(start-1:N-1, :);
Matrix(:, 2) = y(start-2:N-2, :);
Matrix(:, 3) = u(start-1:N-1, :);
Matrix(:, 4) = u(start-2:N-2, :);

est_val = inv(Matrix'*Matrix)*Matrix'*y(start:N,:)
G = [0, 1; est_val(2), est_val(1)];
H = [0;1];
C = [est_val(4), est_val(3)];
D = 0;

sys1 = ss(G, H, C, D, Ts);

y_offset = ave;
u_offset = 2.1;
% values for Kp and Ki
Kp = 0.7;
Ki = 0.025;
sim('lab5')