function num = prox(a, u)

if u >= max(a, 0)
    num = u - a;
elseif u < min(-a, 0)
    num = u + a;
else
    num = 0;
end
      