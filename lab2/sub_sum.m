function C = sub_sum(A,B) 
%this function is used to find the sum of the two subspaces A and B,
%and returns the matrix C 

C = [A B];
C = orth(C);
end 