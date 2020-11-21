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
p2 = [-1 -2 -3 -20];
K1 = place(A,B,p);
K2 = place(A,B,p2);

% define initial conditions
ic = [-0.5, 0, -pi/4, 0];

%set relative tolerances for integration
options = odeset('RelTol', 1e-7,'AbsTol', 1e-7);

%set Tspan
Tspan = linspace(0,10,1e3);

%integrate for both initial conditions
[t,x]=ode45(@cartPendulum,Tspan,ic,options,parameters,K1,A,B);
[t,x2]=ode45(@cartPendulum,Tspan,ic,options,parameters,K2,A,B);

X1_plt = x(:,1); %x1(t): first column of x
X2_plt = x(:,2); %x2(t): third column of x
X3_plt = x(:,3); %x1(t): first column of z
X4_plt = x(:,4); %x2(t): third column of z

X1o_plt = x2(:,1); %x1(t): first column of x
X2o_plt = x2(:,2); %x2(t): third column of x
X3o_plt = x2(:,3); %x1(t): first column of z
X4o_plt = x2(:,4); %x2(t): third column of z


figure('Name', 'K1');
subplot(3,2,1);
plot(t, X1_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
hold on;

subplot(3,2,2);
plot(t, X2_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
hold on;

subplot(3,2,3);
plot(t, X3_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
hold on;

subplot(3,2,4);
plot(t, X4_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
hold on;

subplot(3,2,5);
plot(t, x);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

legend('state 1 of set 1', 'state 2 of set 1', 'state 3 of set 1', 'state 4 of set 1')
hold on;
% ====================================== %
subplot(3,2,1);
plot(t, X1o_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
legend('state 1 of set 1', 'state 1 of set 2')

subplot(3,2,2);
plot(t, X2o_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
legend('state 2 of set 1', 'state 2 of set 2')

subplot(3,2,3);
plot(t, X3o_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
legend('state 3 of set 1', 'state 3 of set 2')

subplot(3,2,4);
plot(t, X4o_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
legend('state 4 of set 1', 'state 4 of set 2')

subplot(3,2,5);
plot(t, x2);
title('x vs time');
xlabel('t');
ylabel('x');

legend('state 1 of set 1', 'state 2 of set 1', 'state 3 of set 1', 'state 4 of set 1', 'state 1 of set 2', 'state 2 of set 2', 'state 3 set 2', 'state 4 of set 2')

