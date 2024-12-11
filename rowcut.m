function D = rowcut(C,i)
% To remove ith row from C (size: N x M) and
% return D (size: N-1 x M)

[m,n] = size(C);
d1 = C(1:i-1,:);
d2 = C(i+1:m,:);
D = [d1; d2];

