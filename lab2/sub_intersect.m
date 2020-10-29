function D = sub_intersect(A,B)

x = [A -B]
n = null(x)

%need to find the set of vectors that satisfy that line x = l*v
%u*w given the nullspace we found
lambda = [n(1);n(2)];
mu = [n(3); n(4)]; 

%use B to verify answer since B*u should equal A*lambda 
v = B*mu
D = A*lambda; 

end 