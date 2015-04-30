function x = indct( v ) 
% INDCT   Compute the inverse non-uniform discrete cosine transform
% 
%   X = INDCT( V ) compute the INDCT with input V. 
% 
tol = 1e5*norm(v, 2)*eps; 
max_iter = 300; 

[x, flag, relres, iter] = cgs(@(x) ndct(x), v, tol, max_iter, @(x) P(x));

end

function y = P( x ) 
% Preconditioner for CGS. 
y = chebfun.idct(flipud(x), 3); 
y(1) = y(1)/2; 
end