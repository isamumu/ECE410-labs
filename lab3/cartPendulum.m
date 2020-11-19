function xdot = cartPendulum(t,x,parameters) %t-time;x-state;parameters-M,g,l
    %extract pendulum parameters
    M = parameters.M; %extract mass
    g = parameters.g; %extract gravitational constant
    l = parameters.l; %extract length of pendulum
    m = parameters.m;

    %extract pendulum states
    x1o = x(1);
    x2o = x(2);
    x3o = x(3);
    x4o = x(4);
    
    xo = [x1o x2o x3o x4o]
    
    %set control signal (u=0 initially)
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

    % ================= 3. Analysing the Cart-Pendulum Model ======================
    % construct controllability matrix
    Qc = ctrb(Alinear,Blinear)
    rankQc = rank(Qc) % notice this equals 4 
    % since n = 4 = rank(Qc), system (A,B) is controllable

    % ===== obtain those gains =====
    % define desired poles
    p = [-1 -2 -3 -4];
    K1 = place(A,B,p);
    
    u = -K1*xo'
    
    %compute xdot
    xdot = [];
    xdot(1) = x2;
    xdot(2) = ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
    xdot(3) = x4;
    xdot(4) = ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)));
    
    xdot = xdot(:);
end