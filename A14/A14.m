clear all;

delay = 2;
st_appServer = 0.03;
st_DBMS = 0.08;
st_Storage = 0.1;
in_vec = [3, 2, 0, 0];
lambda_tot = sum(in_vec);
X = lambda_tot

service_times = [delay, st_appServer, st_Storage, st_DBMS];

p12 = 0.8;
p10 = 0.2;
p24 = 0.5;
p20 = 0.2;
p23 = 0.3;

transition_matrix = [
     0, p12,   0,      0;
     0, 0,   p23,    p24;
     0, 1, 0,      0;
     0, 1, 0,    0;
    ];

reference_load = [1, 0, 0, 0];

in_vec = [3, 2, 0, 0];
load_k = in_vec / X;

visits = load_k * inv(eye(4) - transition_matrix);

D1 = visits(1)*delay;
D2 = visits(2)*st_appServer;
D3 = visits(3)*st_Storage;
D4 = visits(4)*st_DBMS;



R1 = D1;

U2 = lambda_tot*D2;
R2 = D2/(1-U2);

U3 = lambda_tot*D3;
R3 = D3/(1-U3);

U4 = lambda_tot*D4;
R4 = D4/(1-U4);

R_avg = (R1 + R2 + R3 + R4)
N = R_avg*lambda_tot