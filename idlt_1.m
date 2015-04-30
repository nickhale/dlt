function c = idlt_1(f)
%IDLT  Inverse discrete Legendre transform.

% Nick Hale and Alex Townsend, Feb 2015.

% Stage I:
n = size(f, 1);
[~, w] = legpts(n);                           % Legendre weights
wv = bsxfun(@times, w.', f);                  % Scale by weights
c = ndct_tranpose(wv);                        % NDCT transpose

% Stage II:
c = leg2cheb(c, 'transpose');                 % LEG2CHEB transpose
c = bsxfun(@times, (0:n-1).' + .5, c);        % Scaling

end

function c = ndct_tranpose(f)
%NDCT_TRANSPOSE  Evalute transpose of the NDCT operator.

L = 18;                                       % Number of terms in series
[n, m] = size(f); nn = (0:n-1).';             % Transform size
[~, ~, ~, t_leg] = legpts(n);                 % Legendre grid (in theta)
t_leg = t_leg(end:-1:1); f = f(end:-1:1,:);   % Get right-left ordering
t_cheb = ((0:n-1).'+.5)*pi./n;                % Chebyshev nodes (in theta)
dt = (t_leg - t_cheb);                        % (t_leg - t_cheb)
dt(n:-1:ceil((n+1)/2)) = -dt(1:ceil(n/2));    % More accurate near t = 0
dt = repmat(dt, 1, m); nn = repmat(nn, 1, m); % Support matrix input
NN = 1; scl = 1; sgn = -1;                    % Initialise
c = chebfun.dct( f, 2 );                      % ell = 0 term
for l = 1:(L-1)                               % Remaining terms in sum
    f = dt.*f; NN = nn.*NN; scl = scl*l;      % Update terms
    if ( mod(l, 2) )
        c = c + (sgn/scl)*NN.*dst3_shifted_transpose(f); % Sine terms
    else
        c = c + (sgn/scl)*NN.*chebfun.dct( f, 2 );       % Cosine terms
    end
    if ( l/2 == round(l/2) ),
        sgn = -sgn;                           % +--++--.
    end
end
% Ensure real/imag output for real/imag input:
if ( isreal(f) ), c = real(c); elseif ( isreal(1i*f) ), c = imag(c); end

end

function v = dst3_shifted_transpose(c)
% dst3_shifted_transpose is computed by a DST-II.
v = chebfun.dst( c, 2 );
v = [zeros(1, size(c, 2)) ; v(1:end-1,:)];
end
