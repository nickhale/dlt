function c = idlt_quad(v)
%IDLT_QUAD  Quad precision implementation of IDLT written in C.

% Note: This code is intended for testing purposes only. To avoid complilation 
% problems with MATLAB MEX compiler we read and write the inputs and outputs to 
% file rather than pass them directly.

% Write coefficients to file.
fid = fopen('data.bin', 'w+');
fprintf(fid, '%16.16f\n', v);
fclose(fid);

% Write Legendre nodes and weights to file.
n = size(v, 1);
[x, w] = legpts(n);
fid = fopen('data_x.bin', 'w+');
fprintf(fid, '%16.16f\n', x);
fclose(fid);
fid = fopen('data_w.bin', 'w+');
fprintf(fid, '%16.16f\n', w);
fclose(fid);

if ( ~exist('idlt_quad.out', 'file') )
    % Compile if required.
    !gcc idlt_quad.c -o idlt_quad.out -lquadmath
end

% Call C function:
!./idlt_quad.out

% Read in results from file.
fid = fopen('results.bin', 'r');
c = fscanf(fid, '%f');
fclose(fid);

% Remove tmp files.
!rm data.bin
!rm data_x.bin
!rm data_w.bin
!rm results.bin

end
