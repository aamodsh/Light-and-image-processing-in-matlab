function imout = hpf(im1,ringsize)
ringsize = 1; %Ring to be excluded in the fourier domain
Fim1 = fftshift(fft2(im1));
filtfunc = ones(size(Fim1));
[X,Y] = meshgrid(1:size(Fim1,2),1:size(Fim1,1));
X = X- mean(X(:));
Y = Y - mean(Y(:));
filtfunc = double((X.^2 + Y.^2) > ringsize.^2);
imout = abs(ifft2(ifftshift(filtfunc.*Fim1)));