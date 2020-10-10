% ================= 3. Numerical Integration ======================
% define three pendulum parameters
M = 0.2; %kg
l = 0.15; %m
g = 9.81; %m/s^2

%define structure
parameters.M = M;
parameters.l = l;
parameters.g = g;

%set initial conditions (2x1 column vector)
x0_0 = [0; sqrt(g/l)];  
x0_1 = [0;1.99.*sqrt(g/l)];
%set relative tolerances for integration
options = odeset('RelTol', 1e-7,'AbsTol', 1e-7);

%set Tspan
Tspan = linspace(0,10,1e3);

%integrate for both initial conditions
[t,x]=ode45(@pendulum,Tspan,x0_0,options,parameters);
[t,x_o]=ode45(@pendulum,Tspan,x0_1,options,parameters);

%============= plotting =============%

% first initial conditions
x1_plt = x(:,1); %x1(t): first column of x
x2_plt = x(:,2); %x2(t): second column of x

figure('Name', 'x vs time (w/ first IC)');
subplot(2,1,1);
plot(t, x1_plt); %plot x1
title('x1(t) vs t')
xlabel('t')
ylabel('x(t)')
subplot(2,1,2);
plot(t, x2_plt); %plot x2
xlabel('t')
ylabel('x(t)')
title('x2(t) vs t')

figure('Name', 'x2 vs x1 (w/ first IC)');
plot(x1_plt, x2_plt);
xlabel('x1');
ylabel('x2');
title('x2(t) vs x1(t)');

% second initial conditions
x1_plt_o = x_o(:,1); %x1(t): first column of x
x2_plt_o = x_o(:,2); %x2(t): second column of x

figure('Name', 'x vs time (w/ second IC)');
subplot(2,1,1);
plot(t, x1_plt_o); %plot x1
title('x1(t) vs t')
xlabel('t')
ylabel('x(t)')
subplot(2,1,2);
plot(t, x2_plt_o); %plot x2
xlabel('t')
ylabel('x(t)')
title('x2(t) vs t')

figure('Name', 'x2 vs x1 (w/ second IC)');
plot(x1_plt_o, x2_plt_o);
xlabel('x1');
ylabel('x2');
title('x2(t) vs x1(t)');

% ================= 4. Symbolic Linearization ======================
syms x1 x2 t u m l g real; %declare symbolic real variables

xdot = [x2;-1.*(g/l).*sin(x1)-(1/(m.*l)).*cos(x1)*u]; %symbolic 2x1
y = [x1 0]; %output equation

%equilibrium points
xbar = 0;
ubar = 0;

%xdot=f and y=h
x = [x1 x2];
A = jacobian(xdot,x);
B = jacobian(xdot,u);
A2 = subs(A, x1, xbar) %should sub in u, but expression goes to 0 already
B2 = subs(B, x1, xbar)

C = jacobian(y,x)
D = jacobian(y,u)

%test new equilibria
syms theta real;
xbar2 = theta;
ubar2 = -1.*m*g*tan(theta);
A3 = subs(A, {x1, u}, {xbar2, ubar2}); %should sub in u, but expression goes to 0 already
A3 = simplify(A3)
B3 = subs(B, x1, xbar2)
%check the printed output: matches expression in chapter 5 of text

% ================= 5. Symbolic to Numerical Integration ======================
syms z1 z2 real;
z = [z1;z2];
x = [x1;x2];

zdot = A*(z-xbar) + B*(u-ubar);
Xdot = [xdot;zdot];
Xdot = subs(Xdot, {m, g, l, u}, {parameters.M, parameters.g, parameters.l, ubar});

symvar(Xdot) %confirm only x and z variables exist
augmented_pend = matlabFunction(Xdot,'Vars',{t,[x;z]});

%initial conditions
X_init1 = [0; sqrt(parameters.g/parameters.l); 0; sqrt(parameters.g/parameters.l)];  
X_init2 = [0; 1.99.*sqrt(parameters.g/parameters.l); 0; 1.99.*sqrt(parameters.g/parameters.l)];

[t,X_0]=ode45(augmented_pend,Tspan,X_init1, options);
[t,X_1]=ode45(augmented_pend,Tspan,X_init2, options);

X1_plt = X_0(:,1); %x1(t): first column of x
X2_plt = X_0(:,2); %x2(t): third column of x
Z1_plt = X_0(:,3); %x1(t): first column of z
Z2_plt = X_0(:,4); %x2(t): third column of z

X1o_plt = X_1(:,1); %x1(t): first column of x
X2o_plt = X_1(:,2); %x2(t): third column of x
Z1o_plt = X_1(:,3); %x1(t): first column of z
Z2o_plt = X_1(:,4); %x2(t): third column of z

% plots for first IC
figure('Name', 'sec. 5: IC = sqrt(g/l)');
subplot(2,1,1); %SUBPLOT 1
plot(t, X1_plt); %plot x1
hold on;
plot(t, Z1_plt, 'r');
title('overlap x1(t) and z1(t)');
xlabel('t');
ylabel('x1(t) or z1(t)');
legend('x1(t)', 'z1(t)');
subplot(2,1,2); %SUBPLOT 2
plot(t, X2_plt); %plot x1
hold on;
plot(t, Z2_plt, 'r');
title('overlap x2(t) and z2(t)');
xlabel('t');
ylabel('x2(t) or z2(t)');
legend('x2(t)', 'z2(t)');

figure('Name', 'sec. 5: IC = sqrt(g/l)');
subplot(2,1,1); %SUBPLOT 1
plot(X1_plt, X2_plt); %plot x1
title('nonlinear orbit of (X1,X2)');
xlabel('x1');
ylabel('x2');
subplot(2,1,2); %SUBPLOT 2
plot(Z1_plt, Z2_plt); %plot x1
title('linearized orbit of (Z1,Z2)');
xlabel('z1');
ylabel('z2');

%plots for second IC
figure('Name', 'sec. 5: IC = 1.99*sqrt(g/l)');
subplot(2,1,1); %SUBPLOT 1
plot(t, X1o_plt); %plot x1
hold on;
plot(t, Z1o_plt, 'r');
title('overlap x1(t) and z1(t)');
xlabel('t');
ylabel('x1(t) or z1(t)');
legend('x1(t)', 'z1(t)');
subplot(2,1,2); %SUBPLOT 2
plot(t, X2o_plt); %plot x1
hold on;
plot(t, Z2o_plt, 'r');
title('overlap x2(t) and z2(t)');
xlabel('t');
ylabel('x2(t) or z2(t)');
legend('x2(t)', 'z2(t)');

figure('Name', 'sec. 5: IC = 1.99*sqrt(g/l)');
subplot(2,1,1); %SUBPLOT 1
plot(X1o_plt, X2o_plt); %plot x1
title('nonlinear orbit of (X1,X2)');
xlabel('x1');
ylabel('x2');
subplot(2,1,2); %SUBPLOT 2
plot(Z1o_plt, Z2o_plt); %plot x1
title('linearized orbit of (Z1,Z2)');
xlabel('z1');
ylabel('z2');

% ================= 6. LTI Representations ======================
A_d = subs(A2, {m,g,l}, {parameters.M, parameters.g, parameters.l}); %should sub in u, but expression goes to 0 already
evalA = eval(A_d)

B_d = subs(B2, {m,g,l}, {parameters.M, parameters.g, parameters.l}); %should sub in u, but expression goes to 0 already
evalB = eval(B_d)

evalC = eval(C)
evalD = eval(D)

A = double(A_d);
B = double(B_d);
C = double(C);
D = double(D);

sys = ss(A,B,C,D)
G = tf(sys) %TF matches the one in textbook
%print eigenvalues of A on the console with: eig(A)

% ================= 7. Pendulum Stabilization ======================
%define a controller with gain as -40 (20log(30)), 2 break frequencies with
%one pole at -10 and one zero at -1000. 
%=======================
num = [1,10];
den1 = [1,1000];
c_lead1 = tf(num, den1)
c_lead = c_lead1 *(-30);

%returns 1+c*g as a factored transfer function 
z = zpk(1+c_lead*G);

%since the zeroes of the numerator of z are contained in the open left half
%plane, we know that our controller results in an output that is BIBO
%stable

%set initial conditions 
x_0 = [pi/2; 0;0]  
%set relative tolerances for integration
options = odeset('RelTol', 1e-7,'AbsTol', 1e-7);

%set Tspan
Tspan = linspace(0,200,1e3);

%integrate for both initial conditions
[t,x]=ode45(@controlled_pendulum,Tspan,x_0,options,parameters);
x1_plt = x(:,1); %x1(t): first column of x
x2_plt = x(:,2); %x2(t): second column of x
z_plt = x(:,3);

figure('Name', 'part 7');
subplot(2,1,1); %SUBPLOT 1
plot(x1_plt, x2_plt); %plot x1
title('nonlinear orbit of (X1,X2)');
xlabel('x1');
ylabel('x2');

figure('Name', 'time plots');
plot(t, z_plt, 'r');
title('zed');
xlabel('t');
ylabel('z(t)');
