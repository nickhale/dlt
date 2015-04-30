function Tc = evalTc(x, c)
%EVALTC  Evaluate Chebyshev-Vandermonde matrix times a vector.

N = length(x);
if ( N < 1 )
    Tc = 1 + 0*c;
    return
end

Tm1 = 1+0*x; T = x;                 % T_0 and T_1.
Tc = c(1) + x*c(2);
for n = 1:(N-2)                     % Recurrence relation:
    Tp1 = 2*x.*T - Tm1;
    Tc = Tc + c(n+2)*Tp1;
    Tm1 = T; T = Tp1;
end

end