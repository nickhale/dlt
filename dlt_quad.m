function v = dlt_quad(c)
%DLT_QUAD  Quad precision implementation of DLT written in C.

% Note: This code is intended for testing purposes only. To avoid complilation 
% problems with MATLAB MEX compiler we read and write the inputs and outputs to 
% file rather than pass them directly.

% Write coefficients to file.
fid = fopen('data.bin', 'w+');
fprintf(fid, '%16.16f\n', c);
fclose(fid);

% Write Legendre nodes to file.
n = size(c, 1);
[~, ~, ~, t] = legpts(n);
fid = fopen('data_t.bin', 'w+');
fprintf(fid, '%16.16f\n', t);
fclose(fid);

if ( ~exist('dlt_quad.out', 'file') )
    % Compile if required.
    !gcc dlt_quad.c -o dlt_quad.out -lquadmath
end

% Execute C code.
!./dlt_quad.out

% Read in results.
fid = fopen('results.bin', 'r');
v = fscanf(fid, '%f');
fclose(fid);

end
