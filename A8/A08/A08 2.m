

% Transition rates between different states
% --- High ---
l01 = 0.33;
l02 = 0;
l03 = 0.05;

% --- Medium ---
l12 = 0.4;
l13 = 0.05;
l23 = 0.05;

% --- Low ---
m10 = 0.6;
m20 = 0;
m21 = 1;

% --- Down ---
m30 = 6 * 0.6;
m31 = 6 * 0.3;
m32 = 6 * 0.1;

% Initializing transition rate sums for each state
s1 = -(l01 + l02 + l03); 
s2 = -(m10 + l12 + l13);
s3 = -(m20 + m21 + l23);
s4 = -(m30 + m31 + m32);

% Building the transition rate matrix Q
Q = [ s1,   l01,    l02,    l03;  %0
      m10,  s2,    l12,    l13;  %1  
      m20,  m21,    s3,     l23;  %2
      m30,  m31,    m32,    s4];  %3

% Initial probability vector
p0 = [0, 1, 0, 0];

% Solving the differential equations for the Markov chain
[t, Sol] = ode45(@(t,x) Q'*x, [0 8], p0');

% Plotting the probabilities of each state over time
plot(t, Sol, "-");