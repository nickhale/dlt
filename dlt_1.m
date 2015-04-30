function v = dlt_1(c_leg)
%DLT  Discrete Legendre transform. 
%   V = DLT(C) returns a column vector V such that
%       V(k) = C(1)*P_0(x(k)) + C(2)*P_1(x(k)) + ... C(N)*P_{N-1}(x(k)), 
%   where P_j(x) is the degree j Legendre polynomial and x is the vector of
%   Gauss-Legendre nodes (as returned by x = legpts(size(C,1))).
%
% See also IDLT, NDCT.

% Nick Hale & Alex Townsend, Feb 2015.

% Stage 1: Convert from Legendre to Chebyshev
c_cheb = leg2cheb(c_leg);

% Stage 2: Compute non-uniform DCT.
v = ndct_1(c_cheb);

end
