clear all;

fprintf("\n\nAnalysis 13\n\n");

% Scenario 1

% Terminals at Station 0
terminals = 10;
N = 10; % Number of users

% CPU at Station 1
service_time_CPU = 20 / 1000;
P12 = 0.3;
P13 = 0.6;
P10 = 0.1;

% Disk at Station 2
service_time_disk = 10 / 1000;
P21 = 0.85;
P23 = 0.15;

% RAM at Station 3
service_time_RAM = 3 / 1000;
P32 = 0.25;
P31 = 0.75;

service_times = [terminals, service_time_CPU, service_time_disk, service_time_RAM];

transition_matrix = [
     0, 1,   0,      0;
     0, 0,   P12,    P13;
     0, P21, 0,      P23;
     0, P31, P32,    0;
    ];

reference_load = [1, 0, 0, 0];

visits_station1 = reference_load * inv(eye(4) - transition_matrix);

demands_station = service_times .* visits_station1;

fprintf("Scenario 1\n");
fprintf("Visits vector for the stations:\n");
disp(visits_station1);
fprintf("Demands vector for the stations:\n");
disp(demands_station);

% Scenario 2

lambda_in = 0.3;
influx_vector = [0.3, 0, 0];
total_lambda = sum(influx_vector); 

% CPU
service_time_CPU = 20 / 1000;
P12 = 0.3;
P13 = 0.6;
P1_out = 0.1;

% Disk
service_time_disk = 10 / 1000;
P21 = 0.8;
P23 = 0.15;
P2_out = 0.05;

% RAM
service_time_RAM = 3 / 1000;
P31 = 0.75;
P32 = 0.25;

service_times = [service_time_CPU, service_time_disk, service_time_RAM];

transition_matrix = [
    0, P12, P13;
    P21, 0, P23;
    P31, P32, 0;
];

load_at_each_station = influx_vector / total_lambda;

visits_stations2 = load_at_each_station * inv(eye(3) - transition_matrix);
throughputs = total_lambda * visits_stations2;
demands_stations2 = service_times .* visits_stations2;

fprintf("Scenario 2\n");
fprintf("Visits vector for the stations:\n");
disp(visits_stations2);
fprintf("Demands vector for the stations:\n");
disp(demands_stations2);
fprintf("Throughputs vector for the stations:\n");
disp(throughputs);
