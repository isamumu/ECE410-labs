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
L1 = K';
K = place(A', C', eval2)
L2 = K';

% calculate K for state feedback controller
p = [-5.1 -5.2 -5.3 -5.4]
K = place(A, B, p)

% Simulate this system using ic: (−0.5, 0, −π/4, 0) and (0, 0, 0, 0) for x
% TODO: create a function expressing the dynamics: feedback and estimator
ic1 = [-0.5, 0, -pi/4, 0];
ic2 = [0,0,0,0];
