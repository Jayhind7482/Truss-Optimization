function D = veccut(C,i)
% To remove member i from C and return D with size 1 less.
% This function is useful in imposing the boundary conditions

[m,n] = size(C);
if m == 1
   d1 = C(1:i-1);
   d2 = C(i+1:n);
   D = [d1 d2];
end
if n == 1
   d1 = C(1:i-1);
   d2 = C(i+1:m);
   D = [d1; d2];
end

