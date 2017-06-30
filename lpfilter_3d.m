function H = lpfilter_3d(type, M, N, P, D0, n, Du, Dv, Dw)
%LPFILTER Computes frequency domain lowpass filters
%   H = LPFILTER(TYPE, M, N, D0, n) creates the transfer function of
%   a lowpass filter, H, of the specified TYPE and size (M-by-N).  To
%   view the filter as an image or mesh plot, it should be centered
%   using H = fftshift(H).
%
%   Valid values for TYPE, D0, and n are:
%
%   'ideal'    Ideal lowpass filter with cutoff frequency D0.  n need
%              not be supplied.  D0 must be positive
%
%   'btw'      Butterworth lowpass filter of order n, and cutoff D0.
%              The default value for n is 1.0.  D0 must be positive.
%
%   'gaussian' Gaussian lowpass filter with cutoff (standard deviation)
%              D0.  n need not be supplied.  D0 must be positive.
%
%   'rect' Rect lowpass filter with Du and Dv as cutoffs

% Use function dftuv to set up the meshgrid arrays needed for
% computing the required distances.
[U, V, W] = dftuv_3d(M, N, P);

% Compute the distances D(U, V).
D = sqrt(U.^2 + V.^2 + W.^2);

% Begin fiter computations.
switch type
    case 'ideal'
        H = double(D <=D0);
    case 'btw'
        if nargin == 4
            n = 1;
        end
        H = 1./(1 + (D./D0).^(2*n));
    case 'gaussian'
        H = exp(-(D.^2)./(2*(D0^2)));
    case 'gaussian1dY'
        Vfftshift = fftshift(V);
        
        H = exp(-(Vfftshift.^2)./(2*(Dv^2)));
        H = ifftshift(H);
    case 'gaussian1dX'
        Ufftshift = fftshift(U);
        H = exp(-(Ufftshift.^2)./(2*(Du^2)));
        H = ifftshift(H);
    case 'gaussian1dZ'
        Wfftshift = fftshift(W);
        H = exp(-(Wfftshift.^2)./(2*(Dw^2)));
        H = ifftshift(H);
    case 'btw1dY'
        Vfftshift = fftshift(V);
        H = 1./(1 + (Vfftshift./Dv).^(2*n));
        H = ifftshift(H);
    case 'btw1dX'
        Ufftshift = fftshift(U);
        H = 1./(1 + (Ufftshift./Du).^(2*n));
        H = ifftshift(H);
   case 'btw1dZ'
        Wfftshift = fftshift(W);
        H = 1./(1 + (Wfftshift./Dw).^(2*n));
        H = ifftshift(H);
    case 'truncSincZ'


        Wfftshift = fftshift(W);
        H = abs(sinc(Wfftshift/Dw));
        midZ = (size(H,3)+1)/2;
        secondLobeEnd = midZ + 2*Dw;
        negsecondLobeEnd = midZ - 2*Dw;
        H(:,:,secondLobeEnd+1:end) = 0;
        H(:,:,1:negsecondLobeEnd) = 0;
        H = ifftshift(H);

%         idxw = find(abs(W) <= Dw);
%         idx = intersect(idxw, intersect(idxu, idxv));
%         H(idx) = 1;

        
    case 'rect'
%          Ufftshift = fftshift(U);
%          Vfftshift = fftshift(V);
        if nargin == 7
%             H = zeros(M, N);
             H = zeros(N, M);
            idxu = find(abs(U) <= Du);
            idxv = find(abs(V) <= Dv);
            idx = intersect(idxu, idxv);
            H(idx) = 1;
%             H = ifftshift(H);
            
        end
    otherwise
        error('Unknown filter type.')
end