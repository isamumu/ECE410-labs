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


figure('Name', 'section 4 plots');
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

% ====================================== %
subplot(3,2,1);
plot(t, X1o_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
legend('state 1 of set 1', 'state 1 of set 2')
hold on;

subplot(3,2,2);
plot(t, X2o_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
legend('state 2 of set 1', 'state 2 of set 2')
hold on;

subplot(3,2,3);
plot(t, X3o_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
legend('state 3 of set 1', 'state 3 of set 2')
hold on;

subplot(3,2,4);
plot(t, X4o_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
legend('state 4 of set 1', 'state 4 of set 2')
hold on;

subplot(3,2,5);
plot(t, x2);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

legend('state 1 of set 1', 'state 2 of set 1', 'state 3 of set 1', 'state 4 of set 1', 'state 1 of set 2', 'state 2 of set 2', 'state 3 set 2', 'state 4 of set 2')

% ================= 5. Linear Quadratic Optimal Control ======================
% Goal change q1
q2 = 5;
q21 = 1;
q22 = 2000;
q1 = 0.05;
R1 = 0.5;
R2 = 0.005;
R3 = 10;

q11 = 0.1;
q12 = 0.005

Q1 = [q11 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0]
Q2 = [q12 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0]
Q3 = [q1 0 0 0; 0 0 0 0; 0 0 q21 0; 0 0 0 0]
Q4 = [q1 0 0 0; 0 0 0 0; 0 0 q22 0; 0 0 0 0]
Q5 = [q1 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0]


K3 = lqr(A,B,Q1,R1);
K4 = lqr(A,B,Q2,R1);

K5 = lqr(A,B,Q3,R1);
K6 = lqr(A,B,Q4,R1);

K7 = lqr(A,B,Q5,R2);
K8 = lqr(A,B,Q5,R3);

[t,x1]=ode45(@cartPendulum,Tspan,ic,options,parameters,K3,A,B);
[t,x2]=ode45(@cartPendulum,Tspan,ic,options,parameters,K4,A,B);
[t,x3]=ode45(@cartPendulum,Tspan,ic,options,parameters,K5,A,B);
[t,x4]=ode45(@cartPendulum,Tspan,ic,options,parameters,K6,A,B);
[t,x5]=ode45(@cartPendulum,Tspan,ic,options,parameters,K7,A,B);
[t,x6]=ode45(@cartPendulum,Tspan,ic,options,parameters,K8,A,B);

X11_plt = x1(:,1); %x1(t): first column of x
X12_plt = x1(:,2); %x2(t): third column of x
X13_plt = x1(:,3); %x1(t): first column of z
X14_plt = x1(:,4); %x2(t): third column of z

X21_plt = x2(:,1); %x1(t): first column of x
X22_plt = x2(:,2); %x2(t): third column of x
X23_plt = x2(:,3); %x1(t): first column of z
X24_plt = x2(:,4); %x2(t): third column of z

X31_plt = x3(:,1); %x1(t): first column of x
X32_plt = x3(:,2); %x2(t): third column of x
X33_plt = x3(:,3); %x1(t): first column of z
X34_plt = x3(:,4); %x2(t): third column of z

X41_plt = x4(:,1); %x1(t): first column of x
X42_plt = x4(:,2); %x2(t): third column of x
X43_plt = x4(:,3); %x1(t): first column of z
X44_plt = x4(:,4); %x2(t): third column of z

X51_plt = x5(:,1); %x1(t): first column of x
X52_plt = x5(:,2); %x2(t): third column of x
X53_plt = x5(:,3); %x1(t): first column of z
X54_plt = x5(:,4); %x2(t): third column of z

X61_plt = x6(:,1); %x1(t): first column of x
X62_plt = x6(:,2); %x2(t): third column of x
X63_plt = x6(:,3); %x1(t): first column of z
X64_plt = x6(:,4); %x2(t): third column of z

% ========= plots ==========
% ==== vary q1 ==== %
figure('Name', 'varying q1');
subplot(3,2,1);
plot(t, X11_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
hold on;

subplot(3,2,2);
plot(t, X12_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
hold on;

subplot(3,2,3);
plot(t, X13_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
hold on;

subplot(3,2,4);
plot(t, X14_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
hold on;

subplot(3,2,5);
plot(t, x1);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

subplot(3,2,1);
plot(t, X21_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
legend('state 1 of q1 = 0.1', 'state 1 of q1 = 0.005')
hold on;

subplot(3,2,2);
plot(t, X22_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
legend('state 2 of q1 = 0.1', 'state 2 of q1 = 0.005')
hold on;

subplot(3,2,3);
plot(t, X23_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
legend('state 3 of q1 = 0.1', 'state 3 of q1 = 0.005')
hold on;

subplot(3,2,4);
plot(t, X24_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
legend('state 4 of q1 = 0.1', 'state 4 of q1 = 0.005')
hold on;

subplot(3,2,5);
plot(t, x2);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

legend('state 1 of q1 = 0.1', 'state 2 of q1 = 0.1', 'state 3 of q1 = 0.1', 'state 4 of q1 = 0.1', 'state 1 of q1 = 0.005', 'state 2 of q1 = 0.005', 'state 3 of q1 = 0.005', 'state 4 of q1 = 0.005')

% ==== vary q2 ==== %
figure('Name', 'varying q2');
subplot(3,2,1);
plot(t, X31_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
hold on;

subplot(3,2,2);
plot(t, X32_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
hold on;

subplot(3,2,3);
plot(t, X33_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
hold on;

subplot(3,2,4);
plot(t, X34_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
hold on;

subplot(3,2,5);
plot(t, x3);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

subplot(3,2,1);
plot(t, X41_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
legend('state 1 of q2 = 1', 'state 1 of q2 = 2000')
hold on;

subplot(3,2,2);
plot(t, X42_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
legend('state 2 of q2 = 1', 'state 2 of q2 = 2000')
hold on;

subplot(3,2,3);
plot(t, X43_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
legend('state 3 of q2 = 1', 'state 3 of q2 = 2000')
hold on;

subplot(3,2,4);
plot(t, X44_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
legend('state 4 of q2 = 1', 'state 4 of q2 = 2000')
hold on;

subplot(3,2,5);
plot(t, x2);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

legend('state 1 of q2 = 1', 'state 2 of q2 = 1', 'state 3 of q2 = 1', 'state 4 of q2 = 1', 'state 1 of q2 = 2000', 'state 2 of q2 = 2000', 'state 3 of q2 = 2000', 'state 4 of q2 = 2000')

% ==== vary R ==== %
figure('Name', 'varying R');
subplot(3,2,1);
plot(t, X51_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
hold on;

subplot(3,2,2);
plot(t, X52_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
hold on;

subplot(3,2,3);
plot(t, X53_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
hold on;

subplot(3,2,4);
plot(t, X54_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
hold on;

subplot(3,2,5);
plot(t, x5);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

subplot(3,2,1);
plot(t, X61_plt);
title('x1 vs time');
xlabel('t');
ylabel('x1');
legend('state 1 of R = 0.005', 'state 1 of R = 10')
hold on;

subplot(3,2,2);
plot(t, X62_plt);
title('x2 vs time');
xlabel('t');
ylabel('x2');
legend('state 2 of R = 0.005', 'state 2 of R = 10')
hold on;

subplot(3,2,3);
plot(t, X63_plt);
title('x3 vs time');
xlabel('t');
ylabel('x3');
legend('state 3 of R = 0.005', 'state 3 of R = 10')
hold on;

subplot(3,2,4);
plot(t, X64_plt);
title('x4 vs time');
xlabel('t');
ylabel('x4');
legend('state 4 of R = 0.005', 'state 4 of R = 10')
hold on;

subplot(3,2,5);
plot(t, x6);
title('x vs time');
xlabel('t');
ylabel('x');
hold on;

legend('state 1 of R = 0.005', 'state 2 of R = 0.005', 'state 3 of R = 0.005', 'state 4 of R = 0.005', 'state 1 of R = 10', 'state 2 of R = 10', 'state 3 of R = 10', 'state 4 of R = 10')

% ================= 6. Nonlinear Comparison ======================
ic = [-1 0 pi/4 0];
syms K x1 x2 x3 x4 t;
x = [x1 x2 x3 x4];

u = K*x';

M = parameters.M; %extract mass
g = parameters.g; %extract gravitational constant
l = parameters.l; %extract length of pendulum
m = parameters.m;
xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];
    
Xdot = subs(xdot,{K}, {K7});
CLS = matlabFunction(Xdot,'Vars',{t,x});
[t,X]=ode45(CLS,Tspan,ic,options)
