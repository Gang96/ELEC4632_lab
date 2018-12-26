u_offset = 2;
y_offset = 2.2365;
% system using dataset SysIdenData_lab4.mat
est_val = [1.2512, -0.2724, 0.0147, 0.0412];
% recreate of transfer function and state space model
Ts = 0.75;
G_D = [0, 1; est_val(2), est_val(1)];
H_D = [0;1];
C_D = [est_val(4), est_val(3)];
D_D = 0;

G = G_D';
H = C_D';
C = H_D';
D = D_D;

p0 = 0.92;
p1 = 0.92;
L_ndb = acker(G, H, [p0, p1])
sys3 = ss(G-H*L_ndb, H, C-D*L_ndb, D, Ts)
Gcl = dcgain(sys3)