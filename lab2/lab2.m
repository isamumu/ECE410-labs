%---Part 3.1: Basic Operations on Subspaces---

%define subspaces v and w 
span_V = [1 -1 0 1; 1 1 0 0; 3 1 0 1]';
span_W = [1 0 2 1; 1 0 -2 0]';

%find the basis for both v and w respectively by using the built-in matlab
%function orth (output 1)
V = orth(span_V);
W = orth(span_W);

%call function sub_sum to generate the sum of the two defined subspaces
sum_VW = sub_sum(V,W)

%call function sub_intersect to find the intersection for the two defined
%subspaces (output 2 is line 13 and 17)
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

