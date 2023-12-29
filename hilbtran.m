function x = hilbtran(xr,n)

%HILBTRAN  calculation of Hilbert transform.
%   Xi = HILBTRAN(Xr) computes imaginary part Xi that is the Hilbert  
%   transform of real vector Xr.
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%   If Xr is a matrix, then HILBTRAN operates along the columns of Xr.
%
%   HILBTRAN(Xr,N) computes the N-point Hilbert transform.  Xr is padded with 
%   zeros if it has less than N points, and truncated if it has more.  
%
%   For a discrete-time analytic signal X, the last half of fft(X) is zero, 
%   and the first (DC) and center (Nyquist) elements of fft(X) are purely real.
%
%   EXAMPLE:
%          Xr = [1 2 3 4];
%          Xi = hilbtran(Xr)
%          produces Xi=imag(X)=[1 -1 -1 1] 

%%
if nargin < 2, n = []; end
if ~isreal(xr)
  warning('Imaginary part ignored!')
  xr = real(xr);
end

%% Work along the first nonsingleton dimension
[xr,nshifts] = shiftdim(xr);
if isempty(n)
  n = size(xr,1);
end

%% n-point FFT over columns.
x = ifft(xr,n,1); 

%% Window for IFFT (Oppenheim, Schaffer, Discrete-Time Signal Processing...)
h  = zeros(n,~isempty(x)); % nx1 for nonempty. 0x0 for empty.
if n > 0 && 2*fix(n/2) == n
  % even and nonempty
  h(n/2+1) = 0;
  h(1:n/2) = 1;
  h(n/2+2:n) = -1;
elseif n>0
  % odd and nonempty
  h(1:(n+1)/2) = 1;
  h((n+3)/2:n) = -1;
end

%% Windowed IFFT (real part = 0)
x = (fft(x.*h(:,ones(1,size(x,2)))));
% x*10000
x = imag(x);

%% Back to the original shape conversion.
x = shiftdim(x,-nshifts);
