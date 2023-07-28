%%
F = @(mat)fftshift(fft(mat));
iF = @(mat)ifft(ifftshift(mat));
F2 = @(mat)fftshift(fft2(mat));
iF2 = @(mat)ifft2(ifftshift(mat));
iF2real = @(mat)ifft2(ifftshift(mat),'symmetric');
