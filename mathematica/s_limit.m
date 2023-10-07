s[d_] := Sum[Binomial[i+2, 2], {i, 0, d}];
sx[d_] := Sum[(d + 2 - i) * (i + 1), {i, 0, d+1}] + Sum[(d + 1 - i) * (i + 1), {i, 0, d}];

f[d_] := s[d] / (3 * sx[d] - 2 * s[d]);

Limit[f[x], x -> Infinity];
