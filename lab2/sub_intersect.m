function D = sub_intersect(A,B)
%should the output just be the nullspace of the matrices? tbd  
x = [A -B];
D = null(x);

%we probably need to find the set of vectors that satisfy that line x = l*v
%+ u*w given the nullspace we found
end 