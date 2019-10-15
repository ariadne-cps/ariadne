param N := 20;

var x{1..N} := 2;

minimize f:
sum {i in 1..N} i*(x[i]^2 + x[i]^4);

subject to cons1:
sum {i in 1..N} x[i]^2 = 1;


solve;

display x;
