% ================= 3. Analysing the Cart-Pendulum Model ======================
% define structure
parameters.M = 1.0731;
parameters.m = 0.2300;
parameters.l = 0.3302;
parameters.g = 0.8;

% define equilibrium
xbar = [0 0 0 0]';

% define symbolic variables
syms x1 x2 x3 x4 u M m l g;

xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];

x = [x1 x2 x3 x4];
A = jacobian(xdot,x)
B = jacobian(xdot,u)

Asub = subs(A,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});
Bsub = subs(B,{M,m,l,g}, {parameters.M, parameters.m, parameters.l, parameters.g});
Alinear = subs(Asub, {x1,x2,x3,x4}, {0,0,0,0});
Blinear = subs(Bsub, {x1,x2,x3,x4}, {0,0,0,0});

A = double(Alinear);
B = double(Blinear);

% ================= 4. Controllability and Pole Assignment ======================
% construct controllability matrix
Qc = ctrb(Alinear,Blinear)
rankQc = rank(Qc) % notice this equals 4 
% since n = 4 = rank(Qc), system (A,B) is controllable

% ===== obtain those gains =====
% define desired poles
p = [-1 -2 -3 -4];
K1 = place(A,B,p)

% define initial conditions
ic = [-0.5, 0, -pi/4, 0];

%set relative tolerances for integration
options = odeset('RelTol', 1e-7,'AbsTol', 1e-7);

%set Tspan
Tspan = linspace(0,10,1e3);

%integrate for both initial conditions
[t,x]=ode45(@cartPendulum,Tspan,ic ,options,parameters);

