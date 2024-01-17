



clear all;


% Loop per numero n di giocatori
M = 10000;    

partecipants = zeros(2,M);
winners = 0;
i = 1;
T = 0; % Total time 
Ti = 0; % Total time for the player
s = 1; % current state
while i < M 
    % START
    dT = 0;
    %
    if s == 1
    % 1 con P = 0.70 - If C1 then...
    
        if rand() <= 0.7
        % 1 con P = 0.20 - if _caduto_ then...
            %{
            Cade con exponential distribution; lambda = 0.5;
            dT Exponential;
            LOSE - Wait T = 5 minuti per il prossimo giocatore
            T += dT + 5 min;
            %}
            if rand() <= 0.20
                dT = Exponential_pdf(rand(), 0.5);
                Ti = Ti + dT;
                T = T + Ti + 5;
                partecipants(:, i) = [i, Ti];
                i = i + 1; % next player
                Ti = 0; % reset time 
                s = 1;

        % 2 con P = 0.80 - if !_caduto_ then... 
            %{
            Non cade con Erlang distribution; k = 4, lambda = 1.5;
            dT Erlang;
            T += dT
            * Arrivato su C1 ----------
            %}
            else
                dT = Erlang_pdf(4, 1.5);
                Ti = Ti + dT;
                s = 2;
            end    
           
    % 2 con P = 0.30 - If C2 then...
        else
        % 1 con P = 0.70 - if _caduto_ then...
            %{
            Cade con exponential distribution; lambda = 0.25;
            dT Exponential;
            LOSE - Wait T = 5 minuti per il prossimo giocatore
            T += dT + 5 min;
            %}
            if rand() < 0.70
                dT = Exponential_pdf(rand(), 0.25); 
                Ti = Ti + dT;
                T = T + Ti + 5;
                partecipants(:, i) = [i, Ti];
                i = i + 1; % next player
                Ti = 0; % reset time
                s = 1;

        % 2 con P = 0.30 - if !_caduto_ then... 
            % Non cade con Uniform distribution; a = 3, b = 6;
            else
                dT = Uniform_pdf(rand(), 3, 6);
                Ti = Ti + dT;
                s = 3;
            end
            % Arrivato su C2 ----------
            
        end
    %}
        
    % C1
    elseif s == 2
    %
    
    % 1 con P = 0.5 - If Yellow then...
        if rand() <= 0.5
        % 1 con P = 0.75 if _caduto_ then...
            % Cade con exponential distribution; lambda = 0.4;
            % LOSE - Wait T = 5 minuti per il prossimo giocatore
            if rand () <= 0.75
                dT = Exponential_pdf(rand(), 0.4);
                Ti = Ti + dT;
                T = T + Ti + 5;
                partecipants(:, i) = [i, Ti]; % aggiungo il partecipante all'array
                i = i + 1; % next player
                Ti = 0; % reset time
                s = 1;


        % 2 con P = 0.25 if !_caduto_ then...
            % Non cade con Erlang distribution; k = 3, lambda = 2;
            else
                dT = Erlang_pdf(3, 2);
                Ti = Ti + dT;
                s = 3;
            end
            % Arrivato su C2 --------
        
    % 2 con P = 0.5 - If White then...
        else
        % 1 con P = 0.40 - If _caduto_ then...
            % Cade con exponential distribution; lambda = 0.2;
            % LOSE - Wait T = 5 minuti per il prossimo giocatore
            if rand() <= 0.40
                dT = Exponential_pdf(rand(), 0.2);
                Ti = Ti + dT;
                T = T + Ti + 5;
                partecipants(:, i) = [i, Ti]; % aggiungo il partecipante all'array
                i = i + 1; % next player
                Ti = 0; % reset time
                s = 1;


        % 2 con P = 0.60 - If !_caduto_ then...
            % Non cade con exponential distribution; lambda = 0.15;
            else
                dT = Exponential_pdf(rand(), 0.15);
                Ti = Ti + dT;
                s = 3;
            end
            % Arrivato su C2 --------
        end
    %}
    
    % C2
    elseif s == 3
    %
        
    % 1 con P = 0.40 - If _caduto_ then...
        % Cade con Erlang distribution; k = 5, lambda = 4;
        if rand() <= 0.40 
        
            dT = Exponential_pdf(rand(), 4);
            Ti = Ti + dT;
            T = T + Ti + 5;
            partecipants(:,i) = [i, Ti]; % aggiungo il partecipante all'array
            i = i + 1; % next player
            Ti = 0; % reset time
            s = 1;
            

    % 2 con P = 0.60 - If !_caduto_ then...
        % Non cade con Erlang distribution; k = 5, lambda = 4;
        else
            dT = Exponential_pdf(rand(), 4);
            Ti = Ti + dT;
            s = 4;
        end
        % Arrivato su END -----------
    
    %}
    
    % END
    else
    %

    % WIN - Wait T = 5 minuti per il prossimo giocatore
        partecipants(:, i) = [i, Ti];
        T = T + Ti + 5;
        winners = [winners, i];
        i = i + 1;
        Ti = 0;
        s = 1;
        %fprintf("Win");

    %}
    end
    % Fine loop
end    
    
    %{
    
    1. Draw a state machine of the system.
    2. Compute the winning probability.
    3. Compute the average duration of a game, either from the ENTRANCE to the EXIT, or from
    the ENTRANCE to the fall in the LAVA.
    4. Compute the throughput of the system: how many games per hour the room can
    accommodate.
    
    %}

numWin = size(winners ,2);
totT = sum(partecipants(2, :)) ;
avgGameDur = sum(partecipants(2, :)) / M; % in minuti
Throughput =  60 / avgGameDur;


WinRate = numWin / M;
fprintf("%g\n", WinRate);


function F = Exponential_pdf(x, lambda)

    F = -log(x) / lambda;

    %F = lambda * exp(-lambda * x);

end

function F = Erlang_pdf(k, lambda)

    % F = gampdf(x, k, 1 / lambda);

    x = 0;

    for i = 1:k
    
        x = x + Exponential_pdf(rand(), lambda);

    end

    F = x;

end

function F = Uniform_pdf(x, a, b)

    %{
       if(x < a || x > b)

        F = 0;

    else 

        F = 1 / (b - a);

    end 
    %}

    F = a + ((b - a) * x);

end