function P = evalP(x)
%EVALP  Evaluate Legendre-Vandermonde matrix.

N = length(x);
if ( N < 2 )
    P = 1+0*x;
    return
end

P = zeros(N);                               % Vandermonde matrix.
P(:,1:2) = [1+0*x, x];                      % P_0 and P_1.     
for n = 1:(N-2)                             % Recurrence relation:
    P(:,n+2) = (2-1/(n+1))*(P(:,n+1).*x) - (1-1/(n+1))*P(:,n);  
end

end