
% Scenario 1
x = [0:100000]/10;
% Two stages hyper-exponential distribution with:



%% SCENARIO 1

uniform_el = @()5 + 5*rand();
exponential_el = @() -log(rand())/0.1;

K0 = 100;
maxK = 20000;
M = 1000;
DK = 100;
MaxRelErr = 0.04;

gam = 0.95;


serviceRnd = Erlang();   %ERLANG
arrivalRnd = @()HyperExp(); % HYPER



%d_gamma = 1.96;

d_gamma = norminv((1+gam)/2);
erlang_index = 1;

K = K0;

tA = 0;
tC = 0;

U = 0;
U2 = 0;
R = 0;
R_var = 0;
R2_var = 0;
R2 = 0;

X = 0;
X2 = 0;
X_v = [];

J = 0;
J2 = 0;
newIters = K;
R_variance= 0;
index = 1;

while K < maxK
	for i = 1:newIters
		Bi = 0;
		Wi = 0;
        Xi = 0;
        Ji = 0;
		tA0 = tA;
	    TC= 0;
        erlang_index = 1;
        serviceRnd = Erlang();
		for j = 1:M
			a_ji = arrivalRnd();
			s_ji = serviceRnd(erlang_index);
            erlang_index = erlang_index + 1;
			
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
        R_variance = R_variance + (Ri^2/M - Ri/M);
		R2_var = R2_var + R_variance^2;
        R2 = R2 + Ri^2;
		
		Ti = tC - tA0;
		Ui = Bi / Ti;
		U = U + Ui;
		U2 = U2 + Ui^2;

        Xi = M / Ti;
        X = X + Xi;
        X_v(index) = Xi;
        X2 = X2 + Xi^2;

        Ji = Wi/Ti;
        J = Ji + J;
        J2 = J2 + Ji^2;

        index = index + 1;
	end
	
	Rm = R / K;
	Rs = sqrt((R2 - R^2/K)/(K-1));
	CiR = [Rm - d_gamma * Rs / sqrt(K), Rm + d_gamma * Rs / sqrt(K)];
	errR = 2 * d_gamma * Rs / sqrt(K) / Rm;
    
	Um = U / K;
	Us = sqrt((U2 - U^2/K)/(K-1));
	CiU = [Um - d_gamma * Us / sqrt(K), Um + d_gamma * Us / sqrt(K)];
	errU = 2 * d_gamma * Us / sqrt(K) / Um;
    
    R_var_m = R_variance/K;
    R2_var_m = R2_var/K;
    Rvar_var = ( R2_var_m - R_var_m^2)/(K-1);
    Rs = sqrt(Rvar_var);
    
    CiVR = [R_var_m - d_gamma * Rs / sqrt(K), R_var_m + d_gamma * Rs / sqrt(K)];
    down = R_var_m - d_gamma * Rs / sqrt(K);
    up = R_var_m + d_gamma * Rs / sqrt(K);
    errVR = 2 * (up-down)/(up+down);   
    errVR = 0;

    Xm = X / K;
    X_var = sum((Xm - X_v).^2)/(K-1);
    Xs = sqrt(X_var);
    CiX = [Xm - d_gamma * Xs / sqrt(K), Xm + d_gamma * Xs / sqrt(K)];
    down = Xm - d_gamma * Xs / sqrt(K);
    up = Xm + d_gamma * Xs / sqrt(K);

    errX = 2 * (up-down)/(up+down);
    
    Jm = J / K;
   
    Js = sqrt((J2 - J^2/K)/(K-1));
    CiJ = [Jm - d_gamma * Js / sqrt(K), Jm + d_gamma * Js / sqrt(K)];
	errJ = 2 * d_gamma * Js / sqrt(K) / Jm;

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

%% SCENARIO 2

%function results = compute_data(arrival_distr, service_distr)
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
R_var = 0;
R2_var = 0;
R2 = 0;

X = 0;
X2 = 0;
X_v = [];

J = 0;
J2 = 0;
newIters = K;
R_variance= 0;
index = 1;

while K < maxK
	for i = 1:newIters
		Bi = 0;
		Wi = 0;
        Xi = 0;
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

        R_variance = R_variance + ((Ri^2)/M - Ri/M);
		R2_var = R2_var + R_variance^2;
		
		Ti = tC - tA0;
		Ui = Bi / Ti;
		U = U + Ui;
		U2 = U2 + Ui^2;

        Xi = M / Ti;
        X = X + Xi;
        X_v(index) = Xi;
        X2 = X2 + Xi^2;

        Ji = Wi/Ti;
        J = Ji + J;
        J2 = J2 + Ji^2;

        index = index + 1;
	end
	
	Rm = R / K;
	Rs = sqrt((R2 - R^2/K)/(K-1));
	CiR = [Rm - d_gamma * Rs / sqrt(K), Rm + d_gamma * Rs / sqrt(K)];
	errR = 2 * d_gamma * Rs / sqrt(K) / Rm;
    
	Um = U / K;
	Us = sqrt((U2 - U^2/K)/(K-1));
	CiU = [Um - d_gamma * Us / sqrt(K), Um + d_gamma * Us / sqrt(K)];
	errU = 2 * d_gamma * Us / sqrt(K) / Um;
    
    R_var_m = R_variance/K;
    R2_var_m = R2_var/K;
    Rvar_var = ( R2_var_m - R_var_m^2)/(K-1);
    Rs = sqrt(Rvar_var);
    
    CiVR = [R_var_m - d_gamma * Rs / sqrt(K), R_var_m + d_gamma * Rs / sqrt(K)];
    down = R_var_m - d_gamma * Rs / sqrt(K);
    up = R_var_m + d_gamma * Rs / sqrt(K);
    errVR = 2 * (up-down)/(up+down);

    errVR = 0;
    Xm = X / K;
    X_var = sum((Xm - X_v).^2)/(K-1);
    Xs = sqrt(X_var);
    CiX = [Xm - d_gamma * Xs / sqrt(K), Xm + d_gamma * Xs / sqrt(K)];
    down = Xm - d_gamma * Xs / sqrt(K);
    up = Xm + d_gamma * Xs / sqrt(K);

    errX = 2 * (up-down)/(up+down);
    
    Jm = J / K;
   
    Js = sqrt((J2 - J^2/K)/(K-1));
    CiJ = [Jm - d_gamma * Js / sqrt(K), Jm + d_gamma * Js / sqrt(K)];
	errJ = 2 * d_gamma * Js / sqrt(K) / Jm;

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
%end

function F = Erlang()

    K = 10;
    Lambda = 1.5;
    N = 1000;
    % produce samples
    unif_4 = reshape(rand(N * K, 1),[K, N]);
    erlang_distr = -log(prod(unif_4))./Lambda;
    erlang_distr = transpose(erlang_distr);

    F = erlang_distr;

end


function F = HyperExp()
    
    p = [0.1, 0.9] ;
    lambda = [0.02, 0.2];
    cp = cumsum (p);

    Nc = 2;
    
    
        r = rand();
        for i = 1:Nc
            if r < cp(1, i)
                break;
            end
        end
        F = - log(rand()) / lambda(1, i);
   
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