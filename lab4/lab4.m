% ================= 3. Observability and State Estimation ======================
% include everything from lab 3 regarding pendulum linearization

parameters.M = 1.0731;
parameters.m = 0.2300;
parameters.l = 0.3302;
parameters.g = 9.81;

% define equilibrium
xbar = [0 0 0 0]';

% define symbolic variables
syms x1 x2 x3 x4 u M m l g;
x_set = [x1;x2;x3;x4];
x_set2 = [x1; x3];
% define xdot 
xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];

y = [x1; x3];

% NEW: ADD C for Observability
    
% find the jacobians
x = [x1 x2 x3 x4];
A = jacobian(xdot,x)
B = jacobian(xdot,u)
C = jacobian(y, x)

xbar = [0 0 0 0];

% substitute the equilibrium value
A = subs(A, x, xbar) %should sub in u, but expression goes to 0 already
B = subs(B, x, xbar)
C = subs(C, x, xbar)

% substitute the parameter values
Alinear = subs(A,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});
Blinear = subs(B,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});
Clinear = subs(C,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});

Alinear = subs(Alinear, {x1,x2,x3,x4}, {0,0,0,0});
Blinear = subs(Blinear, {x1,x2,x3,x4}, {0,0,0,0});
Clinear = subs(Clinear, {x1,x2,x3,x4}, {0,0,0,0});

% convert fractions to double decimals
A = double(Alinear);
B = double(Blinear);
C = double(Clinear);

% ================= 3.1 Noiseless State Estimation ======================
eval1 = [-10 -11 -12 -13];
eval2 = [-40 -41 -42 -43];

obsv(A,C) % observability matrix rank is 4, therefore observable 

% calculate K matrices, followed by L based on the rule that L = K'
% Remember that Matlab uses the convention that u = −Kx
K = place(A', C', eval1)
L1 = -K';
K = place(A', C', eval2)
L2 = -K';

% calculate K for state feedback controller
p = [-5.1 -5.2 -5.3 -5.4]
K = -place(A, B, p)


% Simulate this system using ic: (−0.5, 0, −π/4, 0) and (0, 0, 0, 0) for x
% TODO: create a function expressing the dynamics: feedback and estimator
x1 = [-0.5;0;-pi/4;0;0;0;0;0];

Tspan = linspace(0,10,1e3);
options = odeset('RelTol', 1e-7,'AbsTol', 1e-7);

% define the linear dynamics
lin = @(t,x) [A*x(1:4) + B*K*x(1:4); (A + L1*C) * x(5:8) + B*K*x(1:4) - L1*C * x(1:4)];
lin2 = @(t,x) [A*x(1:4) + B*K*x(1:4); (A + L2*C) * x(5:8) + B*K*x(1:4) - L2*C * x(1:4)];

[t,x]=ode45(lin,Tspan,x1, options);
[t,x2]=ode45(lin2,Tspan,x1, options);

xdelta = x(:,5:8) - x(:,1:4);
xdelta2 = x2(:,5:8) - x2(:,1:4);

figure('Name', 'decoupled linear states');
subplot(2,2,1);
plot(t, xdelta(:,1));
hold on; 
plot(t, xdelta2(:,1));
title('i = 1');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

subplot(2,2,2);
plot(t, xdelta(:,2));
hold on; 
plot(t, xdelta2(:,2));
title('i = 2');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

subplot(2,2,3);
plot(t, xdelta(:,3));
hold on; 
plot(t, xdelta2(:,3));
title('i = 3');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

subplot(2,2,4);
plot(t, xdelta(:,4));
hold on; 
plot(t, xdelta2(:,4));

title('i = 4');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

syms x1 x2 x3 x4;

M = parameters.M;
m = parameters.m;
l = parameters.l;
g = parameters.g;

xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];
% define the nonlinear dynamics% something wrong with xdot

lin = @(t,x) [x(2); ((-m.*l.*sin(x(3)).*x(4).*x(4))+(m*g*sin(x(3))*cos(x(3))+K*x(1:4)))./(M+m.*sin(x(3)).*sin(x(3)));
                x(4); ((-m*l.*sin(x(3)).*cos(x(3)).*x(4).*x(4))+(M+m).*g.*sin(x(3))+K*x(1:4).*cos(x(3)))./(l.*(M+m.*sin(x(3)).*sin(x(3))));
                    (A + L1*C) * x(5:8) + B*K*x(1:4) - L1*[x(1);x(3)]];
                
lin2 = @(t,x) [x(2); ((-m.*l.*sin(x(3)).*x(4).*x(4))+(m*g*sin(x(3))*cos(x(3))+K*x(1:4)))./(M+m.*sin(x(3)).*sin(x(3)));
                x(4); ((-m*l.*sin(x(3)).*cos(x(3)).*x(4).*x(4))+(M+m).*g.*sin(x(3))+K*x(1:4).*cos(x(3)))./(l.*(M+m.*sin(x(3)).*sin(x(3))));
                    (A + L2*C) * x(5:8) + B*K*x(1:4) - L2*[x(1);x(3)]];


x_i = [-0.5;0;-pi/4;0;0;0;0;0];
[t,x]=ode45(lin,Tspan,x_i, options);
[t,x2]=ode45(lin2,Tspan,x_i, options);

xdelta = x(:,5:8) - x(:,1:4);
xdelta2 = x2(:,5:8) - x2(:,1:4);

figure('Name', 'decoupled nonlinear states');
subplot(2,2,1);
plot(t, xdelta(:,1));
hold on; 
plot(t, xdelta2(:,1));
title('i = 1');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

subplot(2,2,2);
plot(t, xdelta(:,2));
hold on; 
plot(t, xdelta2(:,2));
title('i = 2');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

subplot(2,2,3);
plot(t, xdelta(:,3));
hold on; 
plot(t, xdelta2(:,3));
title('i = 3');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

subplot(2,2,4);
plot(t, xdelta(:,4));
hold on; 
plot(t, xdelta2(:,4));

title('i = 4');
xlabel('t');
ylabel('position');
hold on;
legend('L1', 'L2');

% ================= 3.2 State Estimation with Measurement Noise ======================
% Using the Matlab randn function, create an array of 1000 2 × 1 random vectors with mean and covariance

% Create an anonymous function W(t) interpolating this array over your integration interval.

% Create new functions simulating the state feedback control system plus the observer

% Simulate the system as before, using both L1 and L2 and produce one figure showing the estimation error

% Calculate the mean squared estimation error (MSE) over the final 5 seconds of each simulation.