function dyn = dynamics(t, ic,K,A,B,C,L) %t-time;x-state;parameters-M,g,l
    % TODO: return a vector containing the state feedback and estimator
    %extract pendulum states
    x1 = ic(1)
    x2 = ic(2)
    x3 = ic(3)
    x4 = ic(4)
    
    x = [x1 x2 x3 x4];
    
    y = C*x'
    xdot = (A+B*K)*x';
    xdotHat = (A+L*C)*x' + B*K*x' + L*y;
    
    dyn = [xdot; xdotHat]
    %dyn = dyn(:)
    
end
