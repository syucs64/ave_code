function  x = prox(v1, v2)

cond1 = (v2 >= max(v1, 0));
cond2 = (v2 < min(-v1, 0));
x = (v2 - v1) .* cond1 + (v2 + v1) .* cond2;

