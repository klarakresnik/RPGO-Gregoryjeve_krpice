function [a1, a2] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2)
% a0, a3, b0, b3 tocke v prostoru (1, 3) velikosti

syms k0 h0 k1 h1

% rešimo enačbo a0 == k0 * b0 + h0 * c0 in definiramo števili k0, h0
eqns = a0 == k0 * b0 + h0 * c0;
S = solve(eqns, [k0, h0]);
k0 = S.k0; 
h0 = S.h0;

% rešimo enačbo a3 == k1 * b3 + h1 * c2 in definiramo števili k1, h1
eqns = a3 == k1 * b3 + h1 * c2;
S = solve(eqns, [k1, h1]);
k1 = S.k1; 
h1 = S.h1;

% definiramo b1, b2 z linearno interpolacijo
b1 = 2 / 3 * b0 + 1 / 3 * b3;
b2 = 1 / 3 * b0 + 2 / 3 * b3;

% definiramo a1, a2
a1 = (k1 - k0) * b0 / 3 + k0 * b1 + 2 * h0 * c1 / 3 + h1 * c0 / 3;
a2 = k1 * b2 - (k1 - k0) * b3 / 3 + h0 * c2 / 3 + 2 * h1 * c1 / 3;

end