function v = dct3(c)
%DCT3  Compute the type-III discrete cosine transform (DCT-III).
%   DCT3(C) returns the DCT-III transform of the vector C. If C is a matrix, the
%   transform is applied to the columns.
%
%   DCTs are scaled in many different ways. We have decided to be consistent
%   with Wikipedia: http://en.wikipedia.org/wiki/Discrete_cosine_transform.

% Nick Hale & Alex Townsend, Feb 2015.

[n, m] = size(c);                             % Get sizes.
persistent w                                  % Weights.
if ( size(w, 1) ~= 2*n )
    % Pre-compute the weight vector:
    w = (exp(-1i*(0:2*n-1)*pi/(2*n))/2).';
    w(n+1) = 0;  w(n+2:end) = -w(n+2:end);
end
u = [c ; ones(1, m) ; c(end:-1:2,:)];         % Mirror the values for FFT.
v = fft(repmat(w, 1, m).*u);                  % Apply weights and compute FFT.
v = v(1:n,:);                                 % Extract first n terms.
if ( isreal(c) ), v = real(v); end            % Ensure real --> real.

end