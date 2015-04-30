function v = dlt(c)
%DLT  Discrete Legendre transform. 
%   V = DLT(C) returns a column vector V such that
%       V(k) = C(1)*P_0(x(k)) + C(2)*P_1(x(k)) + ... C(N)*P_{N-1}(x(k)), 
%   where P_j(x) is the degree j Legendre polynomial and x is the vector of
%   Gauss-Legendre nodes (as returned by x = legpts(size(C,1))).
%
% See also IDLT, NDCT.

% Nick Hale & Alex Townsend, Feb 2015.

v = dlt_1(c);

end
