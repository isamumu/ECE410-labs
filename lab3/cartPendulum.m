function xdot = cartPendulum(t,ic,parameters, K) %t-time;x-state;parameters-M,g,l
  
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
    
    xdot = xdot(:)
end