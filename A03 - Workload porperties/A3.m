trace1 = csvread("Trace1.csv");
trace2 = csvread("Trace2.csv");
trace3 = csvread("Trace3.csv");

Ass_3("TRACE1", trace1);
Ass_3("TRACE2", trace2);
Ass_3("TRACE3", trace3);

function [] = Ass_3(titleT, trace)
    N = length(trace(:,1));

    M1 = sum(trace)/N;
    M2 = sum(trace .^2)/N;
    M3 = sum(trace .^3)/N;
    M4 = sum(trace .^4)/N;

    V = sum((trace-M1) .^2)/N;
    C3 = sum((trace-M1) .^3)/N;
    C4 = sum((trace-M1) .^4)/N;

    SM4 = C4/(V^2);
    

    SD = sqrt(V);
    CV = SD/M1;
    SKEWNESS = C3/(SD^3);
    KURTOSIS = SM4 - 3;

    %CDF: sort the set, then apply i/N
    sortedTrace = [sort(trace)];
    i = 1;

    while i < N
        sortedTrace(i, 2) = i/N;
        i = i + 1;
    end
   
    

    %25th 50th 75th percentiles, remember to use a sorted set!
    h25 = (N-1)*(25/100) + 1;
    h50 = (N-1)*(50/100) + 1;
    h75 = (N-1)*(75/100) + 1;

    if floor(h25) == N
        P25 = trace(N,1);
    elseif floor(h25) < N
        P25 = sortedTrace(floor(h25)) + (h25 - floor(h25)) * (sortedTrace(floor(h25)+1) - sortedTrace(floor(h25)));
    end

    if floor(h50) == N
        P50 = trace(N,1);
    elseif floor(h50) < N
        P50 = sortedTrace(floor(h50)) + (h50 - floor(h50)) * (sortedTrace(floor(h50)+1) - sortedTrace(floor(h50)));
    end

    if floor(h75) == N
        P75 = trace(N,1);
    elseif floor(h75) < N
        P75 = sortedTrace(floor(h75)) + (h75 - floor(h75)) * (sortedTrace(floor(h75)+1) - sortedTrace(floor(h75)));
    end
    
    %Pearson Correlation Coefficient: remember to use the not sorted set
    PCC_vector = [zeros(100,2)];

    m = 1;
    while m <= 100
        
        sum1 = 0;
        sum2 = 0;
        i = 1;
    
        while i < N-m
            sum1 = sum1 + (trace(i) - M1)*(trace(m + i) - M1);
            i = i +  1;
        end
        i = 1;
        while i < N
            sum2 = sum2 + (trace(i) - M1)^2;
             i = i +  1;
        end
        numerator = (1/(N-m))*sum1;
        denominator = (1/N)*sum2;
    
        PCC_vector(m,1) = numerator/denominator;
        PCC_vector(m,2) = m;
        m = m +1;
    end
    
    figure;

    subplot(1,2,1);
    plot(sortedTrace(:, 1), sortedTrace(:, 2), ".");
    title('CDF ' + titleT);
    xlabel("Value");
    ylabel("Probability");

    subplot(1,2,2);
    plot(PCC_vector(:,2), PCC_vector(:,1), "-");
    title('Pearson Correlation Coefficient ' + titleT);
    xlabel("m");
    ylabel("Pearson Coefficient");

    fprintf(1, "%s\n", titleT);
    fprintf(1, "Mean: %g\n", M1);
    fprintf(1, "2nd moment: %g\n", M2);
    fprintf(1, "3rd moment: %g\n", M3);
    fprintf(1, "4th moment: %g\n", M4);

    fprintf(1, "Variance: %g\n", V);
    fprintf(1, "3rd centered moment: %g\n", C3);
    fprintf(1, "4th centered moment: %g\n", C4);

    fprintf(1, "Standardized moment: %g\n", SM4);
    fprintf(1, "Skewness: %g\n", SKEWNESS);

    fprintf(1, "Standard deviation: %g\n", SD);
    fprintf(1, "Coefficient of Variation: %g\n", CV);

    fprintf(1, "Excess Kurtosis: %g\n", KURTOSIS);
    
    
    fprintf(1, "First Quartile: %g\n", P25);
    fprintf(1, "Median: %g\n", P50);
    fprintf(1, "Third Quartile: %g\n", P75);


    fprintf(1, "\n#########################\n\n");


end
