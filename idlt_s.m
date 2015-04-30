function c = idlt_s(f)
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

L = 7;                                        % Number of terms in series
[n, m] = size(f); nn = (0:n-1).';             % Transform size
[~, ~, ~, t_leg] = legpts(n);                 % Legendre grid
t_leg = t_leg(end:-1:1); f = f(end:-1:1,:);   % Get right-left ordering
t_star = (nn+.75)*pi./(n+1/2);                % Shifted equigrid
dt = (t_leg - t_star);                        % (t_leg - t_cheb)
dt(n:-1:ceil((n+1)/2)) = -dt(1:ceil(n/2));    % More accurate near t=0
dt = repmat(dt, 1, m); nn = repmat(nn, 1, m); % Support matrix input
NN = 1; scl = 1; sgn = -1;                    % Initialise
c = dct3_scaled_transpose(f);                 % ell = 0 term
for l = 1:L-1                                 % Remaining terms in sum
    f = dt.*f; NN = nn.*NN; scl = scl*l;      % Update terms
    if ( mod(l, 2) )
        c = c + (sgn/scl)*NN.*dst3_shifted_transpose(f); % Sine terms
    else
        c = c + (sgn/scl)*NN.*dct3_scaled_transpose(f);  % Cosine terms
    end
    if ( l/2 == round(l/2) ),
        sgn = -sgn;                           % +--++--.
    end
end
% Ensure real/imag output for real/imag input:
if ( isreal(f) ), c = real(c); elseif ( isreal(1i*f) ), c = imag(c); end

end

function v = dct3_scaled_transpose(c)
[n, m] = size(c);                       % Get sizes.
persistent w                            % Weights.
if ( size(w, 1) ~= n )
    % Pre-compute the weight vector:
    w = (2*n+1)*exp(.5i*(0:(n-1))*pi/(2*n+1)).';
end
C = zeros(2*n+1, m); C(2:2:end,:) = c;  % Pad with zeros.
v = ifft([ C ; C(end:-1:1,:) ]);        % Mirror values and compute FFT.
v = repmat(w, 1, m).*v(1:n,:);          % Truncate and multiply by weights.
end

function v = dst3_shifted_transpose(c)
[n, m] = size(c);                       % Get sizes.
persistent w1 w2
if ( size(w1, 1) ~= n-1 )
    t = 0.5*pi*(1:n-1)'/(2*n+1);
    w1 = .5*sin( t );
    w2 = -.5i*cos( t );
end
z = zeros(1, m);                        % Zero vector.
C = zeros(2*n, m); C(2:2:end,:) = c;    % Pad with zeros.
C_mirror = C(end:-1:1,:);               % Flipped coefficients.
y1 = fft( [ z ; C ; z ; C_mirror ] );   % Mirror values and compute FFT.
v1 = repmat(w1, 1, m).*y1(2:n,:);       % Truncate and multiply by weights.,
y2 = fft( [ z ; C ; z ; -C_mirror ] );  % Mirror values and compute FFT.
idx = (4*n+2):-1:(3*n+4);               % Extract terms.
v2 = repmat(w2, 1, m).*y2(idx,:);       % Multiply by weights.
v = [z ; v2 - v1];                      % Combine.
end


