function D = sub_intersect(A,B)
%should the output just be the nullspace of the matrices? tbd  
x = [A -B]
D = null(x) 
end 