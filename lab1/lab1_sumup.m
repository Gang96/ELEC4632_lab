%**************** Part a ******************
clc;
close all;
clear all;
load('SysIdenData_StudentVersion.mat');
y_actm = LogData.signals(1).values(:,1);
y_act = LogData.signals(1).values(:,2);
u_act = LogData.signals(2).values;
t = LogData.time;
Ts = t(2)-t(1);

%**************** Part b ******************
figure(1)
subplot(2,1,1)
plot(t, y_act, 'b', t, y_actm, 'r');
grid on;
xlim([0, 700]);
ylim([1 4]);
grid on;
legend('Noise-Reduced Output', 'Measured Output');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('\bfActual Output Signal');

subplot(2,1,2)
plot(t, u_act, 'c');
grid on;
xlim([0, 700]);
ylim([1 3]);
grid on;
legend('Actual Input');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');
title('\bfActual Input Signal');

%**************** Part c ******************
i = 1;
ave = 0;
while(u_act(i+1) == u_act(i))
%     ave = mean(y_act(1:i,:));
    i = i+1;
end
ave=mean(y_act(1:i))
y = y_act-ave;
u = u_act-u_act(1);
u_act(1)

figure(2)
subplot(2,1,1)
plot(t, y, 'r');
grid on;
xlim([0, 700]);
ylim([-2 1]);
grid on;
legend('Actual Output');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('\bfActual Offset-Free Output Signal');

subplot(2,1,2)
plot(t, u, 'c');
grid on;
xlim([0, 700]);
ylim([-0.5 0.5]);
grid on;
legend('Actual Input');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');
title('\bfActual Offset-Free Input Signal');

%**************** Part a ******************
%starting from t = 3, N = 452
start = 3;
N = ceil(size(y, 1)/2);

Matrix(:, 1) = y(start-1:N-1, :);
Matrix(:, 2) = y(start-2:N-2, :);
Matrix(:, 3) = u(start-1:N-1, :);
Matrix(:, 4) = u(start-2:N-2, :);

%**************** Part b ******************
est_val = inv(Matrix'*Matrix)*Matrix'*y(start:N,:)
G = [0, 1; est_val(2), est_val(1)];
H = [0;1];
C = [est_val(4), est_val(3)];
D = 0;
%**************** Part c ******************
sys1 = tf([est_val(3), est_val(4)], [1, -est_val(1), -est_val(2)], Ts)
sys2 = ss(G, H, C, D, Ts)

%**************** Part a ******************
SrcDataToSimulink = [t, u];
A = [1, -est_val(1), -est_val(2)];
B = [est_val(3), est_val(4)];
%y_test = filter(B, A, u);
y_test = lsim(sys2, u, t);
disp('MSE for the whole scequence');
MSE = mean((y_test-y).^2)
disp('MSE for the second half scequence');
MSE3 = mean((y_test(N+1:end)-y(N+1:end)).^2)

%**************** Part b ******************
figure(3)
y_sim = filter(B, A, u(N+1:end,:));
y_act = y(N+1:end,:);
t_tr = t(1:N-1,:);
subplot(2,1,1)
plot(t_tr, y_sim,'--', t_tr, y_act);
grid on;
legend('Simulated Output','Actual Output');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('Offset5-Free Model Verification(2^{nd} Half)');
text(10, 0.7, strcat('MSE = ', num2str(MSE3)));

subplot(2,1,2)
plot(t, y_test, '--', t, y);
grid on;
legend('Simulated Output','Actual Output');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('Offset-Free Model Verification(Entire)');
text(10, 0.7, strcat('MSE = ', num2str(MSE)));

%**************** Part c ******************

%**************** Part 1 ******************
% y(k) = -a1y(k-1)+b1u(k-1)
%starting from t = 2, N = 452

start = 10;
Matrix2(:, 1) = y(start-1:N-1, :);
Matrix2(:, 2) = u(start-1:N-1, :);
est_val2 = inv(Matrix2'*Matrix2)*Matrix2'*y(start:N,:);
sys3 = tf([est_val2(2)], [1, -est_val2(1)])
G = [0, 1; 0, est_val2(1)];
H = [0;1];
C = [0, est_val2(2)];
D = 0;
sys4 = ss(G, H, C, D, Ts)
% sys4 = ss(-est_val(1), est_val(2))
A = [1, -est_val2(1)];
B = [est_val2(2)];
% y_test2 = filter(B, A, u);
y_test2 = lsim(sys4, u, t);
MSE2 = mean((y_test2-y).^2)
figure(4)
plot( t, y, 'r', t, y_test2, 'g.', t, y_test,'b--');
grid on;
legend('Actual Output', '1^{st}Order Model Response', '2^{nd}Order Model Response');
xlabel('Time(sec)');
ylabel('Water Level(V)');
title('Comparison of Different Offset-Free Models');
text(10, 0.65, strcat('MSE1 = ', num2str(MSE2)));
text(10, 0.35, strcat('MSE2 = ', num2str(MSE)));

%**************** Bonus question ******************
%Observing Nyquist Sampling principles, sampling frequency is 2 times greater than frequency of sampled signal. 
%for designing purpose, we set fs>=5fh
%The smaller the sampling time is, the closer tha sampled discrete signal to the continous signal