clear all;

% M/M/1/K and M/M/2/K systems
% A web server receives requests arriving according to a Poisson process of rate lambda = 150 req/min, and serves them with an average service time D = 350 ms. 
% The web server can handle a maximum of 32 requests at a time: if the system is full, new arrivals are discarded.
% Determine:
% 1. The utilization of the system
% 2. The loss probability
% 3. The average number of jobs in the system
% 4. The drop rate
% 6. The average response time
% 7. The average time spent in the queue (waiting for service)

% After 1 year, the load has increased to lambda = 250 req/min, making the current solution no longer applicable. The system administrator adds a second web server and a load balancer: 
% requests enqueues at the load balancer, and then are sent to the first available webserver. 
% Considering the communication time between load balancers and servers to be negligible compared to the service times, determine for this new configuration:
% 1. The total and average utilization of the system 
% 2. The loss probability
% 3. The average number of jobs in the system
% 4. The drop rate
% 6. The average response time
% 7. The average time spent in the queue (waiting for service)

disp("M/M/1/K");

lambda = 150/(60); 
D = 0.35;
K = 32;
rho = lambda*D;

num = rho - rho^(K+1);
den = 1 - rho^(K+1);
util_tot = num/den

num = rho^(K) - rho^(K+1);
den = 1 - rho^(K+1);
P_loss = num/den

drop_rate = lambda*P_loss

num = 1 - (K+1)*rho^(K) + K*rho^(K+1);
den = (1-rho)*(1-rho^(K));
avg_response = D*(num/den)

f1 = rho/(1-rho);
num = (K+1)*rho^(K+1);
den = 1 - rho^(K+1);
avg_jobs = f1 - (num/den)

avg_time_queue = avg_response - D

disp("M/M/2/K");

lambda = 250/(60); 
D = 0.35;
K = 32;
c = 2;
rho = (lambda*D)/c;

num = ((c*rho)^c)*(1 - rho^(K-c+1));
den = factorial(c)*(1-rho);
f1 = num/den;

series = 0;

for i=0:(c-1)
    num = (c*rho)^i;
    den = factorial(i);
    series = series + (num/den);
end

p0 = 1/(f1 + series);

series1 = 0;
series2 = 0;

for i=1:c
    f1 = p0/(factorial(i));
    f2 = (lambda/(1/D))^i;
    pi = f1*f2;

    series1 = series1 + i*pi;
end

for i=(c+1):K

    den = factorial(c)*c^(i-c);
    f1 = p0/den;
    f2 = (lambda/(1/D))^i;
    pi = f1*f2;
    series2 = series2 + pi;

end

series2 = c*series2;

util_tot = series1 + series2
util_avg = util_tot/c

den = factorial(c)*c^(K-c);
f1 = p0/den;
f2 = (lambda/(1/D))^K;
pK = f1*f2;
loss_prob = pK

series1 = 0;
series2 = 0;

for i=1:c
    f1 = p0/(factorial(i));
    f2 = (lambda/(1/D))^i;
    pi = f1*f2;

    series1 = series1 + i*pi;
end

for i=(c+1):K

    den = factorial(c)*c^(i-c);
    f1 = p0/den;
    f2 = (lambda/(1/D))^i;
    pi = f1*f2;
    series2 = series2 + i*pi;

end

avg_jobs = series1 + series2

drop_rate = lambda*pK

den = lambda*(1-pK);
avg_resp = avg_jobs/den

avg_time_queue = avg_resp - D

