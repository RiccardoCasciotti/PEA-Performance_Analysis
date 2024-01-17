
% Consider a server that executes jobs arriving according to a Poisson process of rate lambda = 40 job/s, and serves them with an average service time D = 16 ms.
% Determine:
% 1. The utilization of the system
% 2. The probability of having exactly one job in the system
% 3. The probability of having less than 10 jobs in the system
% 4. The average queue length (jobs not in service)
% 5. The average response time
% 6. The probability that the response time is greater than 0.5 s. 
% 7. The 90 percentile of the response time distribution

% After 1 year, the load has increased to lambda = 90 job/s, making the current solution no longer applicable. The system administrator adds a second server and 
% % load balancer: jobs enqueues at the load balancer, and then are sent to the first available server. 
% Considering the communication time between load balancers and servers to be negligible compared to the service times, determine for this new configuration:
% 1. The total and average utilization of the system
% 2. The probability of having exactly one job in the system 
% 3. The probability of having less than 10 jobs in the system 
% 4. The average queue length (jobs not in service)
% 5. The average response time

clear all;
%% Scenario 1
disp("Scenario 1");
D1 = 0.016;
lambda1 = 40;
rho1 = lambda1*D1;
R1 = 1/((1/D1) - lambda1);
util1 = rho1

avg_resp = 1/( (1/D1) - lambda1 )
P1J = (1-(lambda1/(1/D1)))*(lambda1/(1/D1))
Pless10J = 1- rho1^(10) % or 1- rho1^(11), didn't really understand if it should be less or equal or strictly less
avg_queue =  (rho1^2)/(1 - rho1)
P_resptime = exp(-0.5/R1)
percentile_90 = -log(1-(90/100))*R1


%% Scenario 2
disp("Scenario 2");
lambda2 = 90;

rho2 = lambda2*D1/2;
util_tot = lambda2*D1
util_avg = rho2

p0 = (1-rho2)/(1+rho2);
P1J = 2*p0*rho2

avg_resp = D1/(1-rho2^2)

avg_jsystem = 2*rho2/(1-rho2^2);
avg_jqueue = avg_jsystem - util_tot

Pless10J = p0;

for i = 1:9
    Pless10J = Pless10J + 2*p0*rho2^i;
end

Pless10J






