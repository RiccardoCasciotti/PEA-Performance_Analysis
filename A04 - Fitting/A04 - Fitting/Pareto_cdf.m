function F = Pareto_cdf(x, p)

alpha = p(1);
m = p(2);

F = max(0, min(1, 1 - (m./x).^alpha));

end

