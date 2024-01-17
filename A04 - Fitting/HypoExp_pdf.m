function F = HypoExp_pdf(x, p)
	l1 = p(1);
	l2 = p(2);
	
	F = (x>0) .* (l1*l2/(l1-l2) * (exp(-l2*x) - exp(-l1*x)));
end

%When calling the hypoExponential pdf in the Maximum Likelihood Estimation
%Method (MLE) the starting points shoudl be the one defined in the slides:
%[1/(0.3* M1 ), 1/( 0.7*M1 )]