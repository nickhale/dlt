function T = evalT(x)
%EVALT  Evaluate Chebyshev-Vandermonde matrix.

N = length(x);
if ( N < 2 )
    T = 1+0*x;
    return
end

T = zeros(N);                       % Vandermonde matrix. 
T(:,1:2) = [1+0*x, x];              % T_0 and T_1.
for n = 1:(N-2)                     % Recurrence relation:
    T(:,n+2) = 2*(T(:,n+1).*x) - T(:,n);
end

end