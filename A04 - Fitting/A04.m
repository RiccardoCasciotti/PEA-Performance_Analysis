clear all;

file_one = 'Trace1.csv';
file_two = 'Trace2.csv';


F1 = csvread(file_one);
F2 = csvread(file_two);

size1 = size(F1, 1);
size2 = size(F2, 1);

F1 = sort(F1);
F2 = sort(F2);

%%%%%%%%%%%% FIRST TRACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean
fM = mean(F1);
% Second Moment
sM = mean(F1.^2);
% Third Moment
tM = mean(F1.^3);
% Dev std
sigma = (sM - fM.^2);
% Variance
var = sM - fM.^2;
% Coeficient of Variation
Cv = sqrt(var) / fM;

% Uniform distr
a = 12 * sigma;
a = sqrt(a);
b = fM + (0.5 * a);
a = fM - (0.5 * a);

% Expo distr
lambda = 1 / fM;

% Erlang
Elambda = fM / sigma;
Ek = Elambda.^2 * sigma;


% Weibull
equations = @(x) [x(1) * gamma(1 + 1/x(2)) - fM; x(1)^2 * gamma(1 + 2/x(2)) - sM];
% Stima iniziale per lambda e k
initial_guess = [1, 1];
% Risolvo il sistema di equazioni
solution = fsolve(equations, initial_guess);
Wk = solution(2);
Wl = solution(1);

% Pareto
one = -(fM.^2) + sM;
two = (2 * fM.^2) - (2*sM);
three = - fM.^2;

alpha = ( -two + sqrt(two.^2 - (4*one*three)) ) / (2 * one); 
m = fM * (alpha - 1) / alpha;

% Hyper Exp

p = 0.4;
lambda1 = 0.8 / fM;
lambda2 = 1.2 / fM;

starter = [lambda1, lambda2, p];

Hyperpdf = @(x, l1, l2, p)HyperExp_pdf(x,[l1, l2, p]);

pars1 = mle(F1,"pdf",Hyperpdf,'Start',starter);

% Hypo Exp

lambda1 = 1 / (0.3 * fM);
lambda2 = 1 / (0.7 * fM);

starter = [lambda1,lambda2];

Hypopdf = @(x,l1,l2) HypoExp_pdf(x, [l1,l2]);

pars2 = mle(F1,"pdf",Hypopdf,'Start',starter);

% Plotting


x = [0:60];

% plot(F1, [0:size1-1]/size1,'-'); % Samples
% hold on;
% plot(x, Unif_cdf(x,[a,b]), '-'); % Uniform
% hold on;
% plot(x, Exp_cdf(x,lambda), '-'); % Exponential
% hold on;
% plot(x, Erlang_cdf(x,[Elambda,Ek]), '-'); % Erlang
% hold on;
% plot(x, Weibull_cdf(x, [Wk, Wl]), '-'); % Weibull
% hold on;
% plot(x, Pareto_cdf(x, [alpha, m]), '-'); % Pareto
% hold on;
% plot(x, HyperExp_cdf(x, pars1)); % HyperExponential
% hold on;
% plot(x, HypoExp_cdf(x,pars2)); % HypoExponential
% legend({'Samples','Uniform', 'Exponential', 'Erlang', 'Weibull', 'Pareto', 'HyperExponential', 'HypoExponential'},'Location','southeast')


%%%%%%%%%%%% SECOND TRACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Mean
fM = mean(F2)
% Second Moment
sM = mean(F2.^2)
% Third Moment
tM = mean(F2.^3)
% Dev std
sigma = (sM - fM.^2);
% Variance
var = sM - fM.^2;
% Coeficient of Variation
Cv = sqrt(var) / fM;

% Uniform distr
a = 12 * sigma;
a = sqrt(a);
b = fM + (0.5 * a)
a = fM - (0.5 * a)

% Expo distr
lambda = 1 / fM

% Erlang
Elambda = fM / sigma
Ek = Elambda.^2 * sigma


% Weibull
equations = @(x) [x(1) * gamma(1 + 1/x(2)) - fM; x(1)^2 * gamma(1 + 2/x(2)) - sM];
% Stima iniziale per lambda e k
initial_guess = [1, 1];
% Risolvo il sistema di equazioni
solution = fsolve(equations, initial_guess);
Wk = solution(2)
Wl = solution(1)

% Pareto
one = -(fM.^2) + sM;
two = (2 * fM.^2) - (2*sM);
three = - fM.^2;

alpha = ( -two + sqrt(two.^2 - (4*one*three)) ) / (2 * one)
m = fM * (alpha - 1) / alpha

% Hyper Exp

p = 0.4;
lambda1 = 0.8 / fM;
lambda2 = 1.2 / fM;

starter = [lambda1, lambda2, p];

Hyperpdf = @(x, l1, l2, p)HyperExp_pdf(x,[l1, l2, p]);

pars1 = mle(F2,"pdf",Hyperpdf,'Start',starter);

% Hypo Exp

lambda1 = 1 / (0.3 * fM);
lambda2 = 1 / (0.7 * fM);

starter = [lambda1,lambda2];

Hypopdf = @(x,l1,l2) HypoExp_pdf(x, [l1,l2]);

pars2 = mle(F2,"pdf",Hypopdf,'Start',starter);
  


% Plotting
x = [0:600];

plot(F2, [0:size1-1]/size1,'-'); % Samples
hold on;
plot(x, Unif_cdf(x,[a,b]), '-'); % Uniform
hold on;
plot(x, Exp_cdf(x,lambda), '-'); % Exponential
hold on;
plot(x, Erlang_cdf(x,[Elambda,Ek]), '-'); % Erlang
hold on;
plot(x, Weibull_cdf(x, [Wk, Wl]), '-'); % Weibull
hold on;
plot(x, Pareto_cdf(x, [alpha, m]), '-'); % Pareto
hold on;
plot(x, HyperExp_cdf(x, pars1)); % HyperExponential
%hold on;
%plot(x, HypoExp_cdf(x,pars2)); % HypoExponential


%}


fprintf(1, 'Hypo-exponential distribution is NOT available for Trace2 (CoV < 1)\nHyper-exponential distribution is NOT available for Trace1 (CoV > 1)');











