function [U, V, W] = dftuv_3d(M, N, P)
%DFTUV Computes meshgrid frequency matrices.
%   [U, V] = DFTUV(M, N) computes meshgrid frequency matrices U and
%   V. U and V are useful for computing frequency-domain filter 
%   functions that can be used with DFTFILT.  U and V are both M-by-N.

% Set up range of variables.
u = 0:(M-1);
v = 0:(N-1);
w = 0:(P-1);
% Compute the indices for use in meshgrid
idx = find(u > M/2);
u(idx) = u(idx) - M;
idy = find(v > N/2);
v(idy) = v(idy) - N;
idz = find(w > P/2);
w(idz) = w(idz) - P;
% Compute the meshgrid arrays
[U, V, W] = meshgrid(u, v, w);