function xdot8 = controlled_pendulum(t,x,parameters) 
    M = parameters.M; %mass 
    g = parameters.g; %gravitational constant 
    l = parameters.l; %length of mass 
    
    num = [1,10];
    den1 = [1,1000];
    c_lead1 = tf(num, den1);
    c_lead = c_lead1 *(-30);
    [F,G,H,L] = ssdata(c_lead); % all are 1x1s
    
    x1 = x(1)
    x2 = x(2)
    z = x(3)
    
    y = x1;
   
    %compute zdot as a vector
    zdot = F*z - G*y;
    u = H*z - L*y;
    
    
    %compute xdot as a vector   
    xdot1 = x2;
    xdot2= (-g/l)*sin(x1) - (1/M*l)*cos(x1)*u;
    xdot8 = [xdot1; xdot2; zdot];
end
