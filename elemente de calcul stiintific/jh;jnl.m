function [Q R b] = givens(A, b)
2 [m n] = size(A);
3 Q = eye(m);
4
5 for k = 1 : min(m-1, n)
6 for l = k + 1 : m
7 r = sqrt(A(k,k)ˆ2 + A(l,k)ˆ2);
8 c = A(k,k)/r;
9 s = -A(l,k)/r;
10
11 aux = [c -s; s c]*[A(k,k:n); A(l,k:n)];
12 A(k,k:n) = aux(1, :);
13 A(l,k:n) = aux(2, :);
14
15 aux = [c -s; s c]*[b(k); b(l)];
16 b(k) = aux(1);
17 b(l) = aux(2);
18
19 aux = [c -s; s c]*[Q(k,1:m); Q(l,1:m)];
20 Q(k,1:m) = aux(1, :);
21 Q(l,1:m) = aux(2, :);
22 end
23 end
 Q = Q’;
26 R = A;
27 end
