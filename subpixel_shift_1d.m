function im_mat_shift = subpixel_shift_1d(input_matrix,pixelshift)

%% Fourier domain shifting
%function im_mat_shift = subpixel_shift(input_matrix,pixelshift_y,pixelshift_x)

%Normalize
% im_mat1 = input_matrix/max(max(input_matrix));
im_mat1 = input_matrix;
%Caluculate fft
im_mat_fft = (fft(im_mat1));
% subplot(1,2,1);imagesc(abs(im_mat1).^2);
% subplot(1,2,2);imagesc(abs(im_mat_fft).^2);colorbar;

%Define fourier domain axes
%Corrected sampling
N = length(im_mat1);
center_pixel = N/2+1;
fx = ((1:N)-center_pixel)/N;


%Subpixel shift
p=exp(-1i*2*pi*(fx*pixelshift));
im_mat_shift = ifft((im_mat_fft).*ifftshift(p));