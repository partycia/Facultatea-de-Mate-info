function [valp vecp] = PutereInv(d, s, h, y, maxIter, tol)
  n = length(d);
  vecp = y;
  valp = 0;
  for k = 1 : maxIter
    cond(1) = d(1) * vecp(1) + s(1) * vecp(2);
    cond(n) = s(n - 1) * vecp(n - 1) + d(n) * vecp(n);
    for i = 2 : n - 1
      cond(i) = s(i - 1) * vecp(i - 1) + d(i) * vecp(i) + s(i) * vecp(i + 1);
    end
    if norm(cond - valp * vecp) < tol
      return;
    end
    z = Thomas(s, d - h, s, vecp);
    vecp = z / norm(z);
    s1 = d(n) * vecp(n) ^ 2;
    s2 = 0;
    for i = 1 : n - 1
      s1 = s1 + d(i) * vecp(i) ^ 2;
      s2 = s2 + 2 * s(i) * vecp(i) * vecp(i + 1);
    end
    valp = s1 + s2;
    h = valp;
  end
end