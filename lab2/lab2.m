%---Part 3.1: Basic Operations on Subspaces---

%define subspaces v and w 
span_V = [1 -1 0 1; 1 1 0 0; 3 1 0 1]';
span_W = [1 0 2 1; 1 0 -2 0]';

%(output 1): find the basis for both v and w respectively by using the 
%built-in matlab function orth 

V = orth(span_V);
W = orth(span_W);

%output 2: print sum and intersection of the two subspaces 
%call function sub_sum and sub_intersect to generate the sum and
%intersection of the two defined subspaces

sum_VW = sub_sum(V,W)
int_VW = sub_intersect(V,W)

%---Part 3.2: Linear Transformations and Linear Equations--- 

%define the vectors in the basis 
x1 = [1;1];
x2 = [2;1];

%the change of coordinates is represented by z = P^-1x 
%define P and x 
x = [2;1];
P = [x1, x2];

p_inv = inv(P);
z = p_inv * x; 

%numerically verify that x can be written as an l.c of x1 and x2 and z
%(output 3)
x_verified = z(1)*x1 + z(2)*x2; 

%output 4: find A_hat matrix in terms of A,P,Q 
%define the matrix A and the basis P and Q 
A = [1 0 1 0 0; 2 1 0 -1 0; 0 -1 3 2 2; -1 -2 4 3 2]';
P = [1 1 1 1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Q = [1 1 0 0 1; 1 -1 0 0 0; 0 0 1 1 0; 0 0 1 -1 0; 0 0 0 0 1];

%find an expression for A_hat 
A_prime = A*P; 
Q_inv = inv(Q); 
A_hat = Q_inv*(A_prime)

%output 5: rank, nullity, injectivity/surjectivity 
%get the rank of A matrix 
rank_A = rank(A); 

%get the number of columns and rows in A where m is the # of rows and n is
%the number of columns 
[m, n] = size(A); 

%By the rank nullity theorem: rank(A) + nullity(A) = n 
%and using the fact that nullity(A) = dim(Ker(A)) 
nullity_A = n - rank(A);   

%determine classification of A-matrix 
if (rank_A == m) && (m == n)
    disp('Transformation is an isomorphism')
elseif (rank_A == m)
    disp('Transformation is surjective')
elseif (rank_A == n)
    disp('Transformation is injective') 
else 
    disp('Transformation is neither')
end 

%output 6


%---Part 4: A-Invariance and Representation Thm---- 



