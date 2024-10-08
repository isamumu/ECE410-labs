% ================= 3. Analysing the Cart-Pendulum Model ======================
% define structure
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

% ================= 4. Controllability and Pole Assignment ======================
% construct controllability matrix
Qc = ctrb(Alinear,Blinear)
rankQc = rank(Qc) % notice this equals 4 
% since n = 4 = rank(Qc), system (A,B) is controllable

% ===== obtain those gains =====
% define desired poles
p = [-1 -2 -3 -4];
p2 = [-1 -2 -3 -20];

% calculate the gain values
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

% calculate the feedback expressions
u1 = -K1*x';
u2 = -K2*x2';

% define states to plot
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
title('position vs time');
xlabel('t');
ylabel('position');
hold on;

subplot(3,2,2);
plot(t, X2_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,3);
plot(t, X3_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
hold on;

subplot(3,2,4);
plot(t, X4_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,5);
plot(t, u1);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

% ====================================== %
subplot(3,2,1);
plot(t, X1o_plt);
title('position vs time');
xlabel('t');
ylabel('position');
legend('K1', 'K2')
hold on;

subplot(3,2,2);
plot(t, X2o_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
legend('K1', 'K2')
hold on;

subplot(3,2,3);
plot(t, X3o_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
legend('K1', 'K2')
hold on;

subplot(3,2,4);
plot(t, X4o_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
legend('K1', 'K2')
hold on;

subplot(3,2,5);
plot(t, u2);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

legend('K1', 'K2')
% ================= 5. Linear Quadratic Optimal Control ======================
% Goal change q1 
% define the new parameters
q2 = 5;
q21 = 1;
q22 = 2000;
q1 = 0.05;
R1 = 0.5;
R2 = 0.005;
R3 = 10;

q11 = 0.1;
q12 = 0.005

% define the Q matrices
Q1 = [q11 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0]
Q2 = [q12 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0]
Q3 = [q1 0 0 0; 0 0 0 0; 0 0 q21 0; 0 0 0 0]
Q4 = [q1 0 0 0; 0 0 0 0; 0 0 q22 0; 0 0 0 0]
Q5 = [q1 0 0 0; 0 0 0 0; 0 0 q2 0; 0 0 0 0]

% find the new gain values
K3 = lqr(A,B,Q1,R1);
K4 = lqr(A,B,Q2,R1);

K5 = lqr(A,B,Q3,R1);
K6 = lqr(A,B,Q4,R1);

K7 = lqr(A,B,Q5,R2);
K8 = lqr(A,B,Q5,R3);

% continue with integration for each K value
[t,x1]=ode45(@cartPendulum,Tspan,ic,options,parameters,K3,A,B);
[t,x2]=ode45(@cartPendulum,Tspan,ic,options,parameters,K4,A,B);
[t,x3]=ode45(@cartPendulum,Tspan,ic,options,parameters,K5,A,B);
[t,x4]=ode45(@cartPendulum,Tspan,ic,options,parameters,K6,A,B);
[t,x5]=ode45(@cartPendulum,Tspan,ic,options,parameters,K7,A,B);
[t,x6]=ode45(@cartPendulum,Tspan,ic,options,parameters,K8,A,B);

% obtain the new feedback controllers
u3 = -K3*x1';
u4 = -K3*x2';
u5 = -K3*x3';
u6 = -K3*x4';
u7 = -K3*x5';
u8 = -K3*x6';

% obtain new states for plotting
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
title('position vs time');
xlabel('t');
ylabel('position');
hold on;

subplot(3,2,2);
plot(t, X12_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,3);
plot(t, X13_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
hold on;

subplot(3,2,4);
plot(t, X14_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,5);
plot(t, u3);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

subplot(3,2,1);
plot(t, X21_plt);
title('position vs time');
xlabel('t');
ylabel('position');
legend('q1 = 0.1', ' q1 = 0.005')
hold on;

subplot(3,2,2);
plot(t, X22_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
legend('q1 = 0.1', 'q1 = 0.005')
hold on;

subplot(3,2,3);
plot(t, X23_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
legend('q1 = 0.1', 'q1 = 0.005')
hold on;

subplot(3,2,4);
plot(t, X24_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
legend('q1 = 0.1', 'q1 = 0.005')
hold on;

subplot(3,2,5);
plot(t, u4);
title('feedback vs time');
xlabel('t');
ylabel('u');

legend('q1 = 0.1', 'q1 = 0.005')
hold on;

% ==== vary q2 ==== %
figure('Name', 'varying q2');
subplot(3,2,1);
plot(t, X31_plt);
title('position vs time');
xlabel('t');
ylabel('position');
hold on;

subplot(3,2,2);
plot(t, X32_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,3);
plot(t, X33_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
hold on;

subplot(3,2,4);
plot(t, X34_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,5);
plot(t, u5);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

subplot(3,2,1);
plot(t, X41_plt);
title('position vs time');
xlabel('t');
ylabel('position');
legend('q2 = 1', 'q2 = 2000')
hold on;

subplot(3,2,2);
plot(t, X42_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
legend('q2 = 1', 'q2 = 2000')
hold on;

subplot(3,2,3);
plot(t, X43_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
legend('q2 = 1', 'q2 = 2000')
hold on;

subplot(3,2,4);
plot(t, X44_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
legend('q2 = 1', 'q2 = 2000')
hold on;

subplot(3,2,5);
plot(t, u6);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

legend('q2 = 1', 'q2 = 2000')

% ==== vary R ==== %
figure('Name', 'varying R');
subplot(3,2,1);
plot(t, X51_plt);
title('position vs time');
xlabel('t');
ylabel('position');
hold on;

subplot(3,2,2);
plot(t, X52_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,3);
plot(t, X53_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
hold on;

subplot(3,2,4);
plot(t, X54_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,5);
plot(t, u7);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

subplot(3,2,1);
plot(t, X61_plt);
title('position vs time');
xlabel('t');
ylabel('position');
legend('R = 0.005', 'R = 10')
hold on;

subplot(3,2,2);
plot(t, X62_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
legend('R = 0.005', 'R = 10')
hold on;

subplot(3,2,3);
plot(t, X63_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
legend('R = 0.005', 'R = 10')
hold on;

subplot(3,2,4);
plot(t, X64_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
legend('R = 0.005', 'R = 10')
hold on;

subplot(3,2,5);
plot(t, u8);
title('feedback vs time');
xlabel('t');
ylabel('u');
hold on;

legend('R = 0.005', 'R = 10')

% ================= 6. Nonlinear Comparison ======================
% define the initial conditions
ic = [-1; 0; pi/4; 0]; 
syms K x1 x2 x3 x4 t;
x = [x1; x2; x3; x4];

M = parameters.M; %extract mass
g = parameters.g; %extract gravitational constant
l = parameters.l; %extract length of pendulum
m = parameters.m;

% define the nonlinear system
xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];

% define the feedback controller with particular gain
usub = -K7*x;

% substitute the new feedback controller
Xdot = subs(xdot, u, -K7*x);
% create function to represent the nonlinear system
non_lin = matlabFunction(Xdot,'Vars',{t,x});

% integrate over the nonlinear system
[t,X] = ode45(non_lin,Tspan,ic,options);

% obtain the states for nonlinear system
Xnl1_plt = X(:,1); %x1(t): first column of x
Xnl2_plt = X(:,2); %x2(t): third column of x
Xnl3_plt = X(:,3); %x1(t): first column of z
Xnl4_plt = X(:,4); %x2(t): third column of z

%Plotting linearized system: 
Acl = A-B*K7;
sys = @(t,x) Acl*x; 

[t,xlin] = ode45(sys, Tspan, ic, options); 

xlin1 = xlin(:,1); 
xlin2 = xlin(:,2); 
xlin3 = xlin(:,3); 
xlin4 = xlin(:,4); 

figure('Name', 'Linear and Non-linear Comparison');
subplot(3,2,1);
plot(t, Xnl1_plt);
title('position vs time');
xlabel('t');
ylabel('position');
hold on;

subplot(3,2,2);
plot(t, Xnl2_plt);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,3);
plot(t, Xnl3_plt);
title('phase vs time');
xlabel('t');
ylabel('phase');
hold on;

subplot(3,2,4);
plot(t, Xnl4_plt);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;

subplot(3,2,1);
plot(t, xlin1);
title('position vs time');
xlabel('t');
ylabel('position');
legend('Nonlinear', 'Linear')
hold on;

subplot(3,2,2);
plot(t, xlin2);
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
legend('Nonlinear', 'Linear')
hold on;

subplot(3,2,3);
plot(t, xlin3);
title('phase vs time');
xlabel('t');
ylabel('phase');
legend('Nonlinear', 'Linear')
hold on;

subplot(3,2,4);
plot(t, xlin4);
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
legend('Nonlinear', 'Linear')
hold on;

%plot diff initial conditions: 
ic = [-5; 0; pi/4; 0];
[t,x_ic2] = ode45(non_lin,Tspan,ic,options);

xc1_plt = x_ic2(:,1); %x1(t): first column of x
xc2_plt = x_ic2(:,2); %x2(t): third column of x
xc3_plt = x_ic2(:,3); %x1(t): first column of z
xc4_plt = x_ic2(:,4); %x2(t): third column of z

ic = [-10; 0; pi/4; 0];
[t,x_ic3] = ode45(non_lin,Tspan,ic,options);

xy1_plt = x_ic3(:,1); %x1(t): first column of x
xy2_plt = x_ic3(:,2); %x2(t): third column of x
xy3_plt = x_ic3(:,3); %x1(t): first column of z
xy4_plt = x_ic3(:,4); %x2(t): third column of z

figure('Name', 'Non-linear system w. diff IC');
subplot(3,2,1);
plot(t, Xnl1_plt);
hold on; 
plot(t, xc1_plt); 
hold on; 
plot(t, xy1_plt); 
title('position vs time');
xlabel('t');
ylabel('position');
hold on;
legend('y=-1', 'y=-5', 'y=-10');

subplot(3,2,2);
plot(t, Xnl2_plt);
hold on; 
plot(t, xc2_plt); 
hold on; 
plot(t, xy2_plt); 
title('position velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;
legend('y=-1', 'y=-5', 'y=-10');

subplot(3,2,3);
plot(t, Xnl3_plt);
hold on; 
plot(t, xc3_plt); 
hold on; 
plot(t, xy3_plt); 
title('phase vs time');
xlabel('t');
ylabel('phase');
hold on;
legend('y=-1', 'y=-5', 'y=-10');

subplot(3,2,4);
plot(t, Xnl4_plt);
hold on; 
plot(t, xc4_plt); 
hold on; 
plot(t, xy4_plt); 
title('phase velocity vs time');
xlabel('t');
ylabel('velocity');
hold on;
legend('y=-1', 'y=-5', 'y=-10');


