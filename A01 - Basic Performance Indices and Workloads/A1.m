clear all;

filename = 'barrier.log';


fid = fopen(filename, 'r'); 

IN = [];
OUT = [];

while ~feof(fid)
    line = fgetl(fid); 
    date = line(2:20);
years = str2num(date(1:4));
days =str2num(date(6:8));
hours = str2num(date(10:11));
minutes = str2num(date(13:14));
seconds = str2double(date(16:17));
milliseconds = str2num(date(19:19))*100;
type = line(23:26);

date = datetime(years, 1, 1, hours, minutes, seconds, milliseconds, 'Format', 'yyyy-MM-dd HH:mm:ss.S' ) + days  ;    %duration = years(years-2023+1)*years(years-2023+1) * calyears(1) + days(days) + hours(hours) + minutes(minutes) + seconds(seconds)
if type == "_OUT"
    OUT(end + 1) = posixtime(date);
else
    IN(end + 1) = posixtime(date);
end

end
fclose(fid); 

AC = [IN(:), OUT(:)];

nA = size(AC, 1);
nC = size(AC, 1);

T = AC(nC, 2) - AC(1, 1);

%%%     C         A
Rt = AC(:,2) - AC(:,1);
W = sum(Rt);

Lambda = nA / T;
X = nC / T;

R = W / nC;
N = W / T;

P1 = sum(Rt < 30) / nC;
P2 = sum(Rt < 3*60) / nC;

%To calculate inter arrival times I do the subtraction of 2 column vectors,
%the first one being  and the arrivals vectors without the first element
%second one being  the arrivals one ( without the last element ).

IAT = AC(2:end, 1) - AC(1:end-1, 1);
AV_IAT = sum(IAT)/nA;



evs = [AC(:,1), ones(nA, 1), zeros(nA, 4); AC(:,2), -ones(nC,1), zeros(nC, 4)];
evs = sortrows(evs, 1);  % to put the timestamps in chronological order, this allows us to consider the most recent ( max( A(i) or C(i-1) ) ) for the calculation of the service times
evs(:,3) = cumsum(evs(:,2)); % 3rd: we sum them so that an arrived job is represented by 1 and a completed job by -1, so at each timestamp we know how many jobs are in queue
evs(1:end-1, 4) = evs(2:end,1) - evs(1:end-1,1); % 4th: we have the difference between to sequential timestamps, thus service time of a single job
evs(:,5) = (evs(:,3) > 0) .* evs(:,4); % 5th: greater than zero on col 3 allows us to consider only positive queues:
% 1. This guarantees that if we encounter a timestamp of a completed
% service, this service is for sure a previous one
%
% , we then multiply the number of jobs in the queue for col 4
evs(:,6) = evs(:,3) .* evs(:,4); % 6th: same as above considering also negative queues

B = sum(evs(:,5));
U = B / T;
S = B / nC;

service_times = evs(:,4);

i = 1;
tot_0 = 0;
tot_1 = 0;
tot_2 = 0;

while i < (nA*2)-1
    if evs(i, 3) == 0
        tot_0 = tot_0 + evs(i+1,1) - evs(i,1);
    elseif evs(i, 3) == 1
        tot_1 = tot_1 + evs(i+1,1) - evs(i,1);
    elseif evs(i, 3) == 2
        tot_2 = tot_2 + evs(i+1,1) - evs(i,1);
    end
   i = i + 1;
end

P3 = length(service_times(service_times > 60))/nC;
P4 = length(IAT(IAT < 60))/nA;
P5 = tot_0/T;
P6 = tot_1/T;
P7 = tot_2/T;

fprintf(1, "Average Number of jobs: %g\n", N);
fprintf(1, "Utilization: %g\n", U);
fprintf(1, "Average Service Time: %g\n", S);
fprintf(1, "Average Inter-arrival Time: %g\n", AV_IAT);
fprintf(1, "Arrival Rate: %g, Throughput %g\n", Lambda, X);
fprintf(1, "Average Response Time: %g\n", R);
fprintf(1, "Pr(R < 30 secs): %g\n", P1);
fprintf(1, "Pr(R < 3 mins): %g\n", P2);
fprintf(1, "Pr(Service Time > 1 min): %g\n", P3);
fprintf(1, "Pr(Inter-arrival Time < 1 min): %g\n", P4);
fprintf(1, "Pr(0 jobs in queue): %g\n", P5);
fprintf(1, "Pr(1 jobs in queue): %g\n", P6);
fprintf(1, "Pr(2 jobs in queue): %g\n", P7);


