clear all;

delay = 40;
st_appServer = 0.05;
st_Storage = 0.002;
st_DBMS = 0.08;
st_D1 = 0.08;
st_D2 = 0.120;
N = 80;

service_times = [delay, st_appServer, st_Storage, st_DBMS, st_D1, st_D2];


p21 = 0.1;
p23 = 0.4;
p24 = 0.5;

p35 = 0.6;
p36 = 0.4;


transition_matrix = [
     0, 1,   0,      0, 0, 0;
     0, 0,   p23,    p24, 0, 0;
     0, 0, 0,      0, p35, p36;
     0, 1, 0,    0, 0, 0;
     0, 1, 0,    0, 0, 0;
     0, 1, 0,    0, 0, 0;
    ];

reference_load = [0, 1, 0, 0, 0, 0];
visits = reference_load * inv(eye(6) - transition_matrix);
demands = service_times.*visits;

Nk = [0, 0, 0, 0, 0, 0];
Rk = [0, 0, 0, 0, 0, 0];
X = 0;
for n = 1:N 
    Rk = [demands(1), demands(2)*(1+Nk(2)), demands(3)*(1+Nk(3)), demands(4)*(1+Nk(4)), demands(5)*(1+Nk(5)), demands(6)*(1+Nk(6))];
    den = delay + sum(Rk);
    X = n/den;
    Nk = X*Rk;
end

X = N/(sum(Rk)+delay)
R_sys = (N/X) - delay
U = X .* demands