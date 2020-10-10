function xdot = pendulum(t,x,parameters) %t-time;x-state;parameters-M,g,l
    %extract pendulum parameters
    M = parameters.M; %extract mass
    g = parameters.g; %extract gravitational constant
    l = parameters.l; %extract length of pendulum
    %extract pendulum states
    x1 = x(1);
    x2 = x(2);
    %set control signal (u=0 initially)
    u = 0;
    
    %compute xdot
    xdot = [];
    xdot(1) = x2;
    xdot(2) = -1.*(g/l).*sin(x1)-(1/(M.*l)).*cos(x1)*u
    xdot = xdot(:);
end
