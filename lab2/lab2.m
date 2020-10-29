%---Part 3.1: Basic Operations on Subspaces---

%define subspaces v and w 
span_V = [1 -1 0 1; 1 1 0 0; 3 1 0 1]';
span_W = [1 0 2 1; 1 0 -2 0]';

%find the basis for both v and w respectively by using the built-in matlab
%function orth
V = orth(span_V);
W = orth(span_W);

%call function sub_sum to generate the sum of the two defined subspaces
sum_VW = sub_sum(V,W);

%call function sub_intersect to find the intersection for the two defined
%subspaces 
int_VW = sub_intersect(V,W);