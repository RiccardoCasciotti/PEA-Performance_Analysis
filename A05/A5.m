
%To generate the random number we are going to use the Linear Congruential
%Generator, which works by computing the random numbers using the mod
%operation on m and on the random number generated at timestamp i-1.
% We divided by m-1 since we don't need to include 1: [0,1)




m = 2^32;
a = 1664525;
c = 1013904223;
seed = 521191478;


random_seq = zeros(1, 10000);
random_seq(1) = seed;
uniform_distr = zeros(1, 10000);
i = 2;
uniform_distr(1) = random_seq(1)/(m-1);
while i < 10001
    random_seq(i) = LCG(random_seq(i-1));
    uniform_distr(i) = random_seq(i)/(m-1);
    i = i + 1;
end

uniform_distr(1)
uniform_distr(2)
uniform_distr(3)


%N1 = 10000 samples of an Exponential distribution of rate lambda = 0.1

lambda = 0.1;
exponential_distr = -log(uniform_distr)./lambda;

figure;
plot(sort(exponential_distr), [1:10000]/10000,[0:250]/10, Exp_cdf([0:250]/10,[lambda] ), "-");
legend("Empirical", "Formula Generated ");
title("Exponential distribution CDF");
xlim([0,25]);
hold on
%N2 = 10000 samples of a Pareto distribution with parameters a = 1.5, m = 5
a = 1.5;
m = 5;

% for i=1:10000
%     val = m/(uniform_distr(i)^(1/a));
%     if val >= m
%         pareto_distr(i) = val;
%     else
%         pareto_distr(i) = 0;
%     end
% end
pareto_distr = m./(uniform_distr.^(1/a));
figure;
plot(sort(pareto_distr), [1:10000]/10000,[0:250]/10, Pareto_cdf([0:250]/10,[a, m] ), "-");
legend("Empirical", "Formula Generated ");
title("Pareto distribution CDF");
xlim([0,25]);

%N3 = 2500 samples of an Erlang distribution with k = 4, and lambda = 0.4
unif_4 = reshape(uniform_distr,[4,2500]);
lambda = 0.4;
erlang_distr = -log(prod(unif_4))./lambda;
figure;

plot(sort(erlang_distr), [1:2500]/2500,[0:250]/10, Erlang_cdf([0:250]/10,[4, lambda] ), "-");
legend("Empirical", "Formula Generated ");
title("Erlang distribution CDF");
xlim([0,25]);

%N4 = 5000 samples of a Hypo-Exponential distribution with rates lambda1 = 0.5, lambda2 = 0.125
% to calculate I proceed by putting the 10000 samples of the uniform into a
% 2 row vector:
lambda1 = 0.5;
lambda2 = 0.125;
unif_2 = reshape(uniform_distr,[5000,2]);
exp_distr1 = -log(unif_2(:,1))./lambda1;
exp_distr2 = -log(unif_2(:,2))./lambda2;

hypo_exp_distr = exp_distr1 + exp_distr2;

figure;
plot(sort(hypo_exp_distr), [1:5000]/5000,[0:250]/10, HypoExp_cdf([0:250]/10,[lambda1, lambda2] ), "-");
title("Hypo-Exponential distribution CDF");
legend("Empirical", "Formula Generated ");
xlim([0,25]);
hold on

%plot(sort(hypo_exp_distr), [1:5000]/5000, ".");

%N5 = 5000 samples of a Hyper-Exponential distribution with rates  lambda1 = 0.5, lambda2 = 0.05,
%p1 = 0.55



p = [0.55, 0.45] ;
lambda = [0.5, 0.05];
cp = cumsum (p);

Nc = 2;

for j=1:5000
    r = uniform_distr(j);
    for i = 1:Nc
        if r < cp(1, i) 
            break;
        end
    end
    hyper_exp_distr(j,1) = - log(uniform_distr(j+5000)) / lambda(1, i);
end
figure;
plot(sort(hyper_exp_distr), [1:5000]/5000,[0:250]/10, HyperExp_cdf([0:250]/10,[lambda(1), lambda(2), p(1)] ), "-");
title("Hyper-Exponential distribution CDF");
legend("Empirical", "Formula Generated ");
xlim([0,25]);
hold on


%The previous distributions correspond to the length of a file, expressed in GB. The charge for storing and making available each file by a provider is:
% a. 0.01 $/GB if the file is less than 10 GB
% b. 0.02 $/GB if the file is greater than 10 GB


files = [file_cost(exponential_distr), file_cost(pareto_distr), file_cost(erlang_distr), file_cost(hypo_exp_distr), file_cost(hyper_exp_distr)]

%file_cost(exponential_distr)

function result = file_cost(distr)

    stop = length(distr);
    result = 0;
    
    for k=1:stop
        if distr(k) >= 10
            result = result + distr(k)*0.02;
        else
            result = result + distr(k)*0.01;
        end
    end
    
end

function res = LCG(num)
    a = 1664525;
    c = 1013904223;
    m = 2^32;
    res = mod(a*num + c, m);
end

function F = HypoExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	
	F = (x>0) .* min(1,max(0,1 - l2/(l2-l1) * exp(-l1*x) + l1/(l2-l1) * exp(-l2*x)));
end
% input vector x
%p(1), first rate
%p(2), second rate
%p(3), probablity of choosing one of the two rates

function F = HyperExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	p1 = p(3);
	
	F = max(0,1 - p1 * exp(-l1*x) - (1-p1) * exp(-l2*x));
end

%p indicates the rate
%x is the vector of data to be plotted
function F = Exp_cdf(x, p)
	l = p(1); %it indicates the rate
	
	F = max(0,1 - exp(-l*x));
end



function F = Erlang_cdf(x, p)

	k = p(1); 
    l = p(2); 
	F =  gamcdf(x, k, 1/l);
end


