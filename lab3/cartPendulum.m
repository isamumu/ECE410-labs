function xdot = cartPendulum(t,ic,parameters) %t-time;x-state;parameters-M,g,l
    
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
    
    % define desired poles
    p = [-1 -2 -3 -4];
    K1 = place(A,B,p)
    K = double(K1)
    
    %extract pendulum states
    x1 = ic(1);
    x2 = ic(2);
    x3 = ic(3);
    x4 = ic(4);
    
    x = [x1 x2 x3 x4];
    
    %set control signal (u=0 initially)
    u = -K*x';
    
    %extract pendulum parameters
    M = parameters.M; %extract mass
    g = parameters.g; %extract gravitational constant
    l = parameters.l; %extract length of pendulum
    m = parameters.m;
    
    %compute xdot
    xdot = [];
    xdot(1) = x2;
    xdot(2) = ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
    xdot(3) = x4;
    xdot(4) = ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)));
    
    xdot = xdot(:);
end