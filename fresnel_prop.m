function [Ef,x,y,Fx,Fy,H] = fresnel_prop(E0,ps,lam,z)
%function Ef = fresnel_prop(E0,ps,lam,z)
%Function Input: Initial field in x-y, wavelength lam, no of sample points N, pixel size in um say,z value
%Function output: Final field in x-y after Fresnel Propagation    
%(ref pg 67,J Goodman, Introduction to Fourier Optics)

%%
%Spatial Sampling
%xsize = 10^3;ysize= 10^3;           %Grid size (Pixel size should be of the order of the wavelength; the sharp edge adds noise in high frequencies)
[M,N] = size(E0);
xsize =  ps*N; ysize = ps*M;

%N = 1001;                          %No of sample points in space (preferably odd to include zero)
x = linspace(-xsize/2,xsize/2,N);   %N point sampling over xsize
y = linspace(-ysize/2,ysize/2,M);

%Proper way of creating frequency axis
wx =2*pi*(0:(N-1))/N; %Create unshifted default omega axis
fx = 1/ps*(wx-pi*(1-mod(N,2)/N))/2/pi; %Shift zero to centre - for even case, pull back by pi, for odd case by pi(1-1/N)

wy =2*pi*(0:(M-1))/M; %Create unshifted default omega axis
fy = 1/ps*(wy-pi*(1-mod(N,2)/N))/2/pi; %Shift zero to centre - for even case, pull back by pi, for odd case by pi(1-1/N)

[Fx,Fy] = meshgrid(fx,fy);

%%
%Point spread function h=H(kx,ky) 
%H = sqrt(z*lam)*exp(1i*pi*lam*z*(Fx.^2+Fy.^2));
H = exp(1i*2*pi/lam*z)*exp(1i*pi*lam*z*(Fx.^2+Fy.^2));
E0fft = fftshift(fft2(E0));                 %Centred about zero (as fx and fy defined to be centred around zero)
G = H.*E0fft;
g = ifft2(ifftshift(G));                    %Output after deshifting the fourier transform
Ef=g;
