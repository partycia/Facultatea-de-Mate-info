function r = Power(a, pw)
  x = eig(a);
  n = length(x);
  r = zeros(n);
  q = zeros(n);
  for i = 1 : n
    b = a - (x(i) * eye(n));
    res = zeros(3, 1);
    X = Jacobi(b, res, [0;0;0], 0.01, 100);
    q(:, i) = X;
    b;
  end
  p = q^-1;
  d = p * a * q;
  r = q * d^pw * p;
  end