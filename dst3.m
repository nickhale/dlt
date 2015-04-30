function v = dst3(c)
%DST3  Compute the type-III discrete sine transform (DST-III).
%   DST3(C) returns the DST-III transform of the vector C. If C is a matrix, the
%   transform is applied to the columns.
%
%   DSTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_sine_transform.

% Nick Hale & Alex Townsend, Feb 2015.
[n, m] = size(c);                             % Get sizes.
persistent w                                  % Weights.
if ( size(w, 1) ~= 2*n )
    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*n-1)*pi/(2*n))/2).';
    w([1,n+1]) = [1, 0];  
    w(n+2:end) = -w(n+2:end);
end
u = [c(end:-1:1,:) ; ones(1, m) ; c(1:end-1,:)];% Mirror the values for FFT.
if ( m > 1 ), wu = repmat(w, 1, m).*u;        % Apply weights.
else          wu = w.*u; end
v = fft(wu);                                  % Compute weighted FFT.
v = v(1:n,:);                                 % Extract first n terms.
v(2:2:end,:) = -v(2:2:end,:);
if ( isreal(c) ), v = real(v); end            % Ensure real --> real.
end