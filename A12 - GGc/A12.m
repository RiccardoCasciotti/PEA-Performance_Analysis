
clear all;
% First Year:
fprintf("\nFirst Year\n\n");

% Jobs Arrival: Poisson Process
lambda1 = 10; % jobs per second

% Duration of Each Job: Hyper-Exponential Distribution
u1 = 50; % jobs per second
u2 = 5; % jobs per second
p1 = 0.8; 
p2 = 1 - p1;

D = (p1 / u1) + (p2 / u2); % Hyper-exp mean 
secondMoment = 2 * ((p1 / u1^2) + (p2 / u2^2));
var1 = secondMoment - D^2;

% Calculating average waiting time w
w = lambda1 * secondMoment / 2;

rho = lambda1 * D;

U = rho;
fprintf("Utilization of the system: %g\n", U);
R1 = D + (w / (1 - rho));
fprintf("Average response time (exact): %g\n", R1);
N1 = rho + (w * (lambda1 / (1 - rho)));
fprintf("Average number of jobs in the system (exact): %g\n", N1);

% After One Year:
fprintf("\nAfter One Year\n\n");

c = 3;

% Traffic: Erlang Distribution 
k = 5; % stages
lambda2 = 240; % jobs per second

mean = k / lambda2; % mean of Erlang
sigma = sqrt(k) / lambda2;
var2 = k / lambda2^2;

% Coefficient of variation for Erlang
cv = sigma / mean;
cv2 = var2 / mean^2;
% Coefficient of variation for Hyper-Exp
ca = sqrt(var1) / D;
ca2 = ca^2;
fact = (cv2 + ca2) / 2;

T = k / lambda2;
rho1 = D / (3 * T);

% Average utilization of the system
U1 = rho1;
fprintf("Average utilization of the system: %g\n", U1);
% Using the M/M/3 formula (L15)
first = D / (c * (1 - rho1));
second = (1 - rho1) * (factorial(c) / (c * rho1)^c);
k0 = (c * rho1)^0 / factorial(0);
k1 = (c * rho1)^1 / factorial(1);
k2 = (c * rho1)^2 / factorial(2);
sum = k0 + k1 + k2;

% Theta: average time spent in the queue
theta = first / (1 + (second * sum));

% Approximate average response time
R2 = D + (fact * theta);
fprintf("Approximate average response time: %g\n", R2);
% Approximate average number of jobs in the system
N2 = R2 / T;
fprintf("Approximate average number of jobs in the system: %g\n", N2);