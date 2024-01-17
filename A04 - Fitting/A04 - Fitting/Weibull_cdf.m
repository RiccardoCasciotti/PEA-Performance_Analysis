function F = Weibull_cdf(x, p)

k = p(1);
lambda = p(2);

F = 1 - exp(-(x/lambda).^k);
end

