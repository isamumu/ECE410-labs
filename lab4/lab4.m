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

% define xdot 
xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];

% NEW: ADD C for Observability
    
% find the jacobians
x = [x1 x2 x3 x4];
A = jacobian(xdot,x)
B = jacobian(xdot,u)

xbar = [0 0 0 0];

% substitute the equilibrium value
A = subs(A, x, xbar) %should sub in u, but expression goes to 0 already
B = subs(B, x, xbar)

% substitute the parameter values
Alinear = subs(A,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});
Blinear = subs(B,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});

Alinear = subs(Alinear, {x1,x2,x3,x4}, {0,0,0,0});
Blinear = subs(Blinear, {x1,x2,x3,x4}, {0,0,0,0});

% convert fractions to double decimals
A = double(Alinear);
B = double(Blinear);

% ================= 3.1 Noiseless State Estimation ======================
eval1 = [-10 -11 -12 -13];
eval2 = [-40 -41 -42 -43];

