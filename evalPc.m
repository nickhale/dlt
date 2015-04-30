function Pc = evalPc(x, c)
%EVALTC  Evaluate Legendre-Vandermonde matrix times a vector.

N = length(x);
if ( N < 1 )
    Pc = 1 + 0*c;
    return
end

Pm1 = 1+0*x; P = x;                 % P_0 and P_1.
Pc = c(1) + x*c(2);
for n = 1:(N-2)                     % Recurrence relation:
    Pp1 = (2-1/(n+1))*(P.*x) - (1-1/(n+1))*Pm1;
    Pm1 = P; P = Pp1;
    Pc = Pc + c(n+2)*P;
end

end