
% Scenario 1
x = [0:100000]/10;
% Two stages hyper-exponential distribution with:
lambda1 =0.02;
lambda2 =0.2;
p1 =0.1;

p = [lambda1, lambda2, p1];
hyper_distr = HyperExp_cdf(x, p);

%Erlang with:
k=10;
lambda=1.5;
p = [k, lambda];
erlang_distr = Erlang_cdf(x, p);

s1 = compute_data(hyper_distr, erlang_distr);

%Scenario 2

uniform_el = @()5 + 5*rand();
exponential_el = @() -log(rand())/0.1;



function results = compute_data(arrival_distr, service_distr)
%% FUNZIONE PER CALCOLARE I RISULTATI, PRESETTATA SULLO SCENARIO 2 PER PROVARE
K0 = 100;
maxK = 20000;
M = 1000;
DK = 100;
MaxRelErr = 0.04;

gam = 0.95;


serviceRnd = @()5+5*rand();   % TO BE MODIFIED: UNIFORME
arrivalRnd = @() -log(rand())/0.1; % ESPONENZIALE



%d_gamma = 1.96;

d_gamma = norminv((1+gam)/2);

K = K0;

tA = 0;
tC = 0;

U = 0;
U2 = 0;
R = 0;
R2 = 0;

X = 0;
X2 = 0;

J = 0;
J2 = 0;
newIters = K;

while K < maxK
	for i = 1:newIters
		Bi = 0;
		Wi = 0;
        Ji = 0;
		tA0 = tA;
	    TC= 0;
		for j = 1:M
			a_ji = arrivalRnd();
			s_ji = serviceRnd();
			
			tC = max(tA, tC) + s_ji;
            TC = tC + TC;
			ri = tC - tA;
			Rd((i-1)*M+j,1) = ri;
	
			tA = tA + a_ji;
			
			Bi = Bi + s_ji;
			
			Wi = Wi + ri;
		end
		
		Ri = Wi / M;
		R = R + Ri;
		R2 = R2 + Ri^2;
		
		Ti = tC - tA0;
		Ui = Bi / Ti;
		U = U + Ui;
		U2 = U2 + Ui^2;
      
        %%% MIO
        Si = Bi/TC;
        Xi = Ui / Si;
        X = X + Xi;
        X2 = X2 + Xi^2;

        Ji = Wi/Ti;
        J = Ji + J;
        J2 = J2 + Ji^2;
        %%%%%
	end
	
	Rm = R / K;
	Rs = sqrt((R2 - R^2/K)/(K-1));
	CiR = [Rm - d_gamma * Rs / sqrt(K), Rm + d_gamma * Rs / sqrt(K)];
	errR = 2 * d_gamma * Rs / sqrt(K) / Rm;

	Um = U / K;
	Us = sqrt((U2 - U^2/K)/(K-1));
	CiU = [Um - d_gamma * Us / sqrt(K), Um + d_gamma * Us / sqrt(K)];
	errU = 2 * d_gamma * Us / sqrt(K) / Um;
    
    % QUESTO  L'HO FATTO IO
    Rs = sqrt((R - Rm^2)/(K-1));
    Rvm = Rs / M;
    CiVR = [Rs - d_gamma * Rs / sqrt(K), Rs + d_gamma * Rs / sqrt(K)];
    errVR = 2 * d_gamma * Rs / sqrt(K) / Rvm;
    errVR = 0;

    Xm = X / K;
    %%% QUESTA E´ LA FORMULA SULLE SLIDE
    Xs = sqrt((X - Xm^2)/(K-1));
    CiX = [Xm - d_gamma * Xs / sqrt(K), Xm + d_gamma * Xs / sqrt(K)];
    errX = 2 * d_gamma * Xs / sqrt(K) / Xm;
    %%% SETTATO A ZERO PER NON FAR TERMINARE LA COMPUTAZIONE, DA ELIMINARE
    %%% SE VUOI CALCOLARE THROUGHTPUT
    errX = 0;

    Jm = J / K;
    %%% QUESTA E´ LA FORMULA DEL PROF
    Js = sqrt((J2 - J^2/K)/(K-1));
    CiJ = [Jm - d_gamma * Js / sqrt(K), Jm + d_gamma * Js / sqrt(K)];
	errJ = 2 * d_gamma * Js / sqrt(K) / Jm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if errR < MaxRelErr && errU < MaxRelErr && errX < MaxRelErr && errJ < MaxRelErr && errVR < MaxRelErr
		break;
	else
		K = K + DK;
		newIters = DK;
	end
end

if errR < MaxRelErr && errU < MaxRelErr && errX < MaxRelErr && errJ < MaxRelErr && errVR < MaxRelErr 
	fprintf(1, "Maximum Relative Error reached in %d Iterations\n", K);
else
	fprintf(1, "Maximum Relative Error NOT REACHED in %d Iterations\n", K);
end	

fprintf(1, "Utilization in [%g, %g], with %g confidence. Relative Error: %g\n", CiU(1,1), CiU(1,2), gam, errU);
fprintf(1, "Resp. Time in [%g, %g], with %g confidence. Relative Error: %g\n", CiR(1,1), CiR(1,2), gam, errR);
fprintf(1, "Throughput in [%g, %g], with %g confidence. Relative Error: %g\n", CiX(1,1), CiX(1,2), gam, errX);
fprintf(1, "Jobs in [%g, %g], with %g confidence. Relative Error: %g\n", CiJ(1,1), CiJ(1,2), gam, errJ);
fprintf(1, "Variance of Response time in [%g, %g], with %g confidence. Relative Error: %g\n", CiVR(1,1), CiVR(1,2), gam, errVR);

results = [];
end

function F = Erlang_cdf(x, p)

	k = p(1); 
    l = p(2); 
	F =  gamcdf(x, k, 1/l);
end

function F = HyperExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	p1 = p(3);
	
	F = max(0,1 - p1 * exp(-l1*x) - (1-p1) * exp(-l2*x));
end

function F = Exp_cdf(x, p)
	l = p(1); %it indicates the rate
	
	F = max(0,1 - exp(-l*x));
end

function F = Unif_cdf(x, p)
	a = p(1);
	b = p(2);
	
	F = max(0, min(1, (x>a) .* (x<b) .* (x - a) / (b - a) + (x >= b)));
end