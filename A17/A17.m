clear all; % Clear workspace

% Arrival rates for different product types
lambdaInA = 2 / 60; % Arrival rate for product type A
lambdaInB = 3 / 60; % Arrival rate for product type B
lambdaInC = 2.5 / 60; % Arrival rate for product type C

% Arrival rates assigned to variables
XA = lambdaInA;
XB = lambdaInB;
XC = lambdaInC;

% Production quantities
S_1A = 10;
S_1B = 4;
S_1C = 6;

% Packaging quantities
S_2A = 12;
S_2B = 3;
S_2C = 6;

% Total arrival rate
X = XA + XB + XC;

% Demands for production and packaging
D1A = S_1A;
D1B = S_1B;
D1C = S_1C;

D2A = S_2A;
D2B = S_2B;
D2C = S_2C;

% Utilization of production and packaging stations
U1A = XA * S_1A;
U1B = XB * S_1B;
U1C = XC * S_1C;

U2A = XA * S_2A;
U2B = XB * S_2B;
U2C = XC * S_2C;

% Total utilization of the two stations
U1 = U1A + U1B + U1C
U2 = U2A + U2B + U2C

% Average number of jobs in the system for each product type
N1A = (XA * D1A) / (1 - U1);
N2A = (XA * D2A) / (1 - U2);
NA = N1A + N2A

N1B = (XB * D1B) / (1 - U1);
N2B = (XB * D2B) / (1 - U2);
NB = N1B + N2B

N1C = (XC * D1C) / (1 - U1);
N2C = (XC * D2C) / (1 - U2);
NC = N1C + N2C

% Average system response time per product type
R1A = N1A / XA;
R2A = N2A / XA;
RA = R1A + R2A

R1B = N1B / XB;
R2B = N2B / XB;
RB = R1B + R2B

R1C = N1C / XC;
R2C = N2C / XC;
RC = R1C + R2C

% Class-independent average system response time
R = (XA * RA / X) + (XB * RB / X) + (XC * RC / X)
