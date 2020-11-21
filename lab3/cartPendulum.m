function xdot = cartPendulum(t,ic,parameters,K,A,B) %t-time;x-state;parameters-M,g,l
  
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
    xdot = A*x' + B*u
    
    xdot = xdot(:)
end