function nonlinearPend = cartPendulum(t,x,parameters,K) %t-time;x-state;parameters-M,g,l
    M = parameters.M; %extract mass
    g = parameters.g; %extract gravitational constant
    l = parameters.l; %extract length of pendulum
    m = parameters.m;
    
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    
    u = -K*x';
    
    xdot = [x2; ((-m.*l.*sin(x3).*x4.*x4)+(m*g*sin(x3)*cos(x3)+u))./(M+m.*sin(x3).*sin(x3));
        x4; ((-m*l.*sin(x3).*cos(x3).*x4.*x4)+(M+m).*g.*sin(x3)+u.*cos(x3))./(l.*(M+m.*sin(x3).*sin(x3)))];
    
    xdot = xdot(:);
    
end