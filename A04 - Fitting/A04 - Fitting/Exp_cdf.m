%p indicates the rate
%x is the vector of data to be plotted
function F = Exp_cdf(x, p)
	l = p(1); %it indicates the rate
	
	F = max(0,1 - exp(-l*x));
end