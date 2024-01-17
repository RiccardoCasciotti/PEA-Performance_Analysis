
%Load the CSV files: the first is the inter-arrival time between two consecutive jobs to a system, 
%                    the second one contains their corresponding service time.
trace1 = csvread("Trace1.csv");
trace2 = csvread("Trace2.csv");
trace3 = csvread("Trace3.csv");


Ass_2("TRACE1", trace1);
Ass_2("TRACE2", trace2);
Ass_2("TRACE3", trace3);


function [] = Ass_2(title, trace)

%Average Response Time: find T by using ai arrival times, then find B(T)
%by using the service times, so we can find U. 
    nC1 = length(trace(:,1));
    
    T1 = sum(trace(:, 1));
    B1 = sum(trace(:, 2));
    U1 = B1/T1;
    
    S1 = B1/nC1;
    
    X1 = nC1/T1;
    
    % We need to calculate the completion time, to do so we can use the formula
    % of the service time si = C(i)^-1 + max(C(i-1)^-1, A(i)^-1). We have the
    % inter-arrival times, so we cumulatively sum them to get the arrival time of the single
    % job.
    
    mat1 = [trace(:,1), trace(:,2), zeros(nC1,3)];
    mat1(:,1) = cumsum(mat1(:,1)); 
    i = 1;
    
    while i < nC1
        if i == 1
            mat1(i, 3) = mat1(i, 2) + mat1(i, 1);
        else
            
            mat1(i, 3) = mat1(i, 2) + max(mat1(i, 1), mat1(i-1, 3));
        end
        i = i + 1;
    end
    
    %Now I calculate the single response times
    
    i = 1;
    idle_count = 0;
    
    while i < nC1
    
        if i < nC1-1 && mat1(i, 3) < mat1(i+1, 1)
            idle_count = idle_count + 1;
        end 
        mat1(i, 4) = mat1(i, 3) - mat1(i, 1);
        i = i + 1;
    end
    
    R1 = sum(mat1(:, 4));
    AV_R1 = R1/nC1;
    
    %When is the system IDLE? When the Arrival of the i+1 job is > than the
    %completion of the i job.
    
    idle_freq = idle_count/T1;
    av_idle = (T1-B1)/idle_count;
    fprintf(1, "%s\n", title);
    fprintf(1, "Average response time: %g\n", AV_R1);
    fprintf(1, "System utilization: %g\n", U1);
    fprintf(1, "Frequency at which the system return idle: %g\n", idle_freq);
    fprintf(1, "Average idle time: %g\n", av_idle);

end

%ASK: to calculate the idle times is it fine if I count the time completion
%time of i i s less than arrival time of i + 1?
