function  x = prox(v1, v2)

x = zeros(size(v1));  
cond1 = (v2 >= max(v1, 0));
x(cond1) = v2(cond1) - v1(cond1);
cond2 = (v2 < min(-v1, 0));
x(cond2) = v2(cond2) + v1(cond2);
