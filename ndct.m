function v = ndct(c)
%NDCT  Evalute a Chebyshev series at Gauss-Legendre nodes.
%   V = NDCT(C) returns a column vector V such that
%       V(k) = C(1)*T_0(x(k)) + C(2)*T_1(x(k)) + ... C(N)*T_{N-1}(x(k)), 
%   where T_j(x) is the degree j Chebyshev polynomial and x is the vector of
%   Gauss-Legendre nodes (as returned by x = legpts(size(C,1))).
%
% See also DLT.

% Nick Hale and Alex Townsend, Feb 2015.

v = ndct_1(c)

end
