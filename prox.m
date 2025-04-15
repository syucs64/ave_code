function  x = prox(v1, v2)
tmp = max(v1, 0);
cond1 = (v2 >= tmp);
cond2 = (v2 < -tmp);
x = (v2 - v1) .* cond1 + (v2 + v1) .* cond2;

