clear all;

M = 1000;
d = 1.96;
a = 0.04; 
K0 = 1;

% First Scenario

%{

% Inter-Arrival Time
% Hyper-exp

Lambda1_H = 0.02;
Lambda2_H = 0.2;
P1_H = 0.1;project_c

% Service Time
% Erlang

K_E = 10;
Lambda_E = 1.5;

% Stopping Criteria

% Utilization
K = K0;
repeat = true;
U = 0;
X = 0;
ART = 0;
AN = 0;
Var = 0;


while repeat
    % Genero i samples
    hyper = Hyper_Samp(M, [Lambda1_H, Lambda2_H, P1_H]);
    erlang = Erlang_Samp(M, [K_E, Lambda_E]);
    % Calcolo il tempo totale nel batch
    Ti = sum(hyper) + erlang(end);
    % Calcolo il busy time nel batch
    Bi = sum(erlang);

    % Utilization
    Ui = Bi / Ti;
    U = [U, Ui];
    xm_1 = sum(U)/K;
    s2_1 = (1 / (K - 1)) * sum((U - xm_1).^2);
    % Trovo il Confidence Interval
    Ul = xm_1 - (d * sqrt(s2_1/K));
    Uu = xm_1 + (d * sqrt(s2_1/K));
    % Calcolo l'errore relativo
    threshold1 = 2 * (Uu-Ul) / (Uu + Ul);

    % Throughput
    Xi = M / Ti;
    X = [X, Xi];
    xm_2 = sum(X)/K;
    s2_2 = (1 / (K - 1)) * sum((X - xm_2).^2);
    % Trovo il Confidence Interval
    Xl = xm_2 - (d * sqrt(s2_2/K));
    Xu = xm_2 + (d * sqrt(s2_2/K));
    % Calcolo l'errore relativo
    threshold2 = 2 * (Xu-Xl) / (Xu + Xl);

    F = [hyper, erlang];

    % Arrival time for i
    evs = [F, zeros(M,4)];
    evs(:,3) = cumsum(evs(:,1));
    evs(1,4) = 0;
    evs(2:end, 4) = evs(1:end-1, 3);
    
    % Completion Time
    mass = zeros(M,1);
    mass(1,1) = evs(1,2);
    for i = 2 : M
        mass(i,1) = max(evs(i,4), mass(i-1,1)) + evs(i, 2);
    end
    evs(:,5) = mass;

    % Response Time
    evs(:,6) = mass(:,1) - evs(:,4);

    % Average Response Time
    ARTi = sum(evs(:,6)) / M;
    ART = [ART, ARTi];
    xm_3 = sum(ART)/K;
    s2_3 = (1 / (K - 1)) * sum((ART - xm_3).^2);
    % Trovo il Confidence Interval
    ARTl = xm_3 - (d * sqrt(s2_3/K));
    ARTu = xm_3 + (d * sqrt(s2_3/K));
    % Calcolo l'errore relativo
    threshold3 = 2 * (ARTu-ARTl) / (ARTu + ARTl);
    
    % Mean
    mean_val = mean(evs(:,6));
    % Second Moment
    s_m = mean(evs(:,6).^2);
    % Variance
    Vari = s_m - mean_val.^2;
    Var = [Var, Vari];
    xm_4 = sum(Var)/K;
    s2_4 = (1 / (K - 1)) * sum((Var - xm_4).^2);
    % Trovo il Confidence Interval
    Varl = xm_4 - (d * sqrt(s2_4/K));
    Varu = xm_4 + (d * sqrt(s2_4/K));
    % Calcolo l'errore relativo
    threshold5 = 2 * (Varu-Varl) / (Varu + Varl);

    ac = [evs(:,4) evs(:,5)]; 

    ev = [ac(:,1), ones(M, 1), zeros(M, 1); ac(:,2), -ones(M,1), zeros(M, 1)];
    ev = sortrows(ev, 1);
    ev(:,3) = cumsum(ev(:,2));
    
    % Average Number of Jobs in the System
    ANi = sum(ev(:,3)) / M;
    AN = [AN, ANi];
    xm_5 = sum(AN) / K;
    s2_5 = (1 / (K - 1)) * sum((AN - xm_5).^2);
    % Trovo il Confidence Interval
    ANl = xm_5 - (d * sqrt(s2_5 / K));
    ANu = xm_5 + (d * sqrt(s2_5 / K));
    % Calcolo l'errore relativo
    threshold4 = 2 * (ANu - ANl) / (ANu + ANl);


    % Incremento il numero di Batch
    K = K + 1;
    % Se il threshold è minore di alpha non devo piu ripetere
    if  (threshold1 <= a) && (threshold2 <= a) && (threshold3 <= a) && (threshold4 <= a) && (threshold5 <= a)
        repeat = false;
    end

end 



%}

% Second Scenario

%

% Exponential

lambdaExp = 0.1;

% Uniform

aU = 5;
bU = 10;

% Stopping Criteria

K = K0;
repeat = true;
U = 0;
X = 0;
ART = 0;
AN = 0;
Var = 0;


while repeat
    % Genero i samples
    exponential = Exponential(M, lambdaExp);
    uniform = Uniform(M, [aU, bU]);
    % Calcolo il tempo totale nel batch
    Ti = sum(exponential) + uniform(end);
    % Calcolo il busy time nel batch
    Bi = sum(uniform);

    % Utilization
    Ui = Bi / Ti;
    U = [U, Ui];
    xm_1 = sum(U)/K;
    s2_1 = (1 / (K - 1)) * sum((U - xm_1).^2);
    % Trovo il Confidence Interval
    Ul = xm_1 - (d * sqrt(s2_1/K));
    Uu = xm_1 + (d * sqrt(s2_1/K));
    % Calcolo l'errore relativo
    threshold1 = 2 * (Uu-Ul) / (Uu + Ul);

    % Throughput
    Xi = M / Ti;
    X = [X, Xi];
    xm_2 = sum(X)/K;
    s2_2 = (1 / (K - 1)) * sum((X - xm_2).^2);
    % Trovo il Confidence Interval
    Xl = xm_2 - (d * sqrt(s2_2/K));
    Xu = xm_2 + (d * sqrt(s2_2/K));
    % Calcolo l'errore relativo
    threshold2 = 2 * (Xu-Xl) / (Xu + Xl);

    F = [exponential, uniform];

    % Arrival time for i
    evs = [F, zeros(M,4)];
    evs(:,3) = cumsum(evs(:,1));
    evs(1,4) = 0;
    evs(2:end, 4) = evs(1:end-1, 3);
    
    % Completion Time
    mass = zeros(M,1);
    mass(1,1) = evs(1,2);
    for i = 2 : M
        mass(i,1) = max(evs(i,4), mass(i-1,1)) + evs(i, 2);
    end
    evs(:,5) = mass;

    % Response Time
    evs(:,6) = mass(:,1) - evs(:,4);

    % Average Response Time
    ARTi = sum(evs(:,6)) / M;
    ART = [ART, ARTi];
    xm_3 = sum(ART)/K;
    s2_3 = (1 / (K - 1)) * sum((ART - xm_3).^2);
    % Trovo il Confidence Interval
    ARTl = xm_3 - (d * sqrt(s2_3/K));
    ARTu = xm_3 + (d * sqrt(s2_3/K));
    % Calcolo l'errore relativo
    threshold3 = 2 * (ARTu-ARTl) / (ARTu + ARTl);
    
    % Mean
    mean_val = mean(evs(:,6));
    % Second Moment
    s_m = mean(evs(:,6).^2);
    % Variance
    Vari = s_m - mean_val.^2;
    Var = [Var, Vari];
    xm_4 = sum(Var)/K;
    s2_4 = (1 / (K - 1)) * sum((Var - xm_4).^2);
    % Trovo il Confidence Interval
    Varl = xm_4 - (d * sqrt(s2_4/K));
    Varu = xm_4 + (d * sqrt(s2_4/K));
    % Calcolo l'errore relativo
    threshold5 = 2 * (Varu-Varl) / (Varu + Varl);

    ac = [evs(:,4) evs(:,5)]; 

    ev = [ac(:,1), ones(M, 1), zeros(M, 1); ac(:,2), -ones(M,1), zeros(M, 1)];
    ev = sortrows(ev, 1);
    ev(:,3) = cumsum(ev(:,2));
    
    % Average Number of Jobs in the System
    ANi = sum(ev(:,3)) / M;
    AN = [AN, ANi];
    xm_5 = sum(AN) / K;
    s2_5 = (1 / (K - 1)) * sum((AN - xm_5).^2);
    % Trovo il Confidence Interval
    ANl = xm_5 - (d * sqrt(s2_5 / K));
    ANu = xm_5 + (d * sqrt(s2_5 / K));
    % Calcolo l'errore relativo
    threshold4 = 2 * (ANu - ANl) / (ANu + ANl);


    % Incremento il numero di Batch
    K = K + 1;
    % Se il threshold è minore di alpha non devo piu ripetere
    if  (threshold1 <= a) && (threshold2 <= a) && (threshold3 <= a) && (threshold4 <= a) && (threshold5 <= a)
        repeat = false;
    end

end 

%}

function F = Hyper_Samp(N, p)

    Lambda1 = p(1);
    Lambda2 = p(2);
    P1 = p(3);

    p = [P1, 1 - P1] ;
    lambda = [Lambda1, Lambda2];
    cp = cumsum (p);
    
    Nc = 2;
    
    distr = rand(2 * N,1);
    
    for j = 1 : N
        r = distr(j);
        for i = 1 : Nc
            if r < cp(1, i) 
                break;
            end
        end
        hyper_exp_distr(j,1) = - log(distr(j + N)) / lambda(1, i);
    end

    F = hyper_exp_distr;

end

function F = Erlang_Samp(N, p)

    K = p(1);
    Lambda = p(2);

    % produce samples
    unif_4 = reshape(rand(N * K, 1),[K, N]);
    erlang_distr = -log(prod(unif_4))./Lambda;
    erlang_distr = transpose(erlang_distr);

    F = erlang_distr;

end

function F = Exponential(N, p)

    lambda = p(1);

    F = -log(rand(N,1))./lambda;


end

function F = Uniform(N, p)

    a = p(1);
    b = p(2);

    unif_distr = a + ((b - a) * rand(N,1));

    F = unif_distr;

end