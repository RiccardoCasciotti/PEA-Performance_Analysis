

clear all;



Q = [(-2*(1/2)*1/2),      (1/2)*1/2,      (1/2)*1/2,      0;
      (1/8)*1/8,         -((1/8)*1/8)-(1/18)*7/8,     0,      (1/18)*7/8;
      (1/3)*1/4,           0,            -((1/3)*1/4)-(1/18)*3/4,      (1/18)*3/4;
      (1/2)*1/2,           0,            (1/2)*1/2,    -(1/2)]
prob_vector = [1/4, 1/4, 1/4, 1/4];
throughput_mat = [0, 1, 1, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 1, 1, 0; ]
%throughput_mat = [1, 0, 0, 1; 1, 0, 0, 1;1, 0, 0, 1;1, 0, 0, 1;];
util = [1, 0, 0, 1]
power = [12, 0.1, 0.1, 12]

% Q = [-3/2,           1/2,      1/2,      0,      0,           0,      1/2,      0;
%       1/3,          -1,        0,    1/3,      0,           0,      0,      1/3;
%       1/8,             0,     -3/8,    1/8,    1/8,           0,      0,      0;
%       0,             1/2,      1/2,   -3/2,      0,         1/2,      0,      0;
%       0,               0,     1/2,       0,   -3/2,         1/2,    1/2,      0;
%       0,               0,       0,    1/18,   1/18,        -1/6,      0,   1/18;
%       1/18,               0,       0,       0,   1/18,           0,   -1/6,  1/18;
%       0,               1/2,       0,       0,      0,         1/2,    1/2,  -3/2];
% throughput_mat = [0, 1, 1, 0, 0, 0, 1, 0;1, 0, 0, 1, 0, 0, 0, 1; 1, 0, 0, 1, 1, 0, 0, 0; 0, 1, 1, 0, 0, 1, 0, 0; 0, 0, 1, 0, 0, 1, 1, 0; 0, 0, 0, 1, 1, 0, 0, 1;1, 0, 0, 0, 1, 0, 0, 1;0, 1, 0, 0, 0, 1, 1, 0;];
% power = [12, 0.1, 0.1, 12, 12, 0.1, 0.1, 12];
% prob_vector = [1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8];
% util = [1, 0, 0, 1, 1, 0, 0, 1 ];


%disp(sum(Q'));
%disp(throughput_mat .* Q);

val = sum((throughput_mat .*Q)');

throughput = (prob_vector*val')*(24*60)

[t, Sol] = ode45(@(t,x) Q'*x, [0 24*60], prob_vector');

%plot(t, Sol, "-");
%legend("S1", "S2", "S3", "S4");



util = [Sol(end,:) * util']
power = [Sol(end,:) * power']
