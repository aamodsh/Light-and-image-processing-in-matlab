function Psi_xy = poisson_solve_symm(func,ps,epsilon)
%function Psi_xy = poisson_solve_symm(func,ps,N,epsilon)
%Solves poisson equation of the form Del2_psi = func, accepts x-y vectors, returns real Psi

[M,N] = size(func);


%Proper way of creating frequency axis
wx =2*pi*(0:(N-1))/N; %Create unshifted default omega axis
fx = 1/ps*(wx-pi*(1-mod(N,2)/N))/2/pi; %Shift zero to centre - for even case, pull back by pi, for odd case by pi(1-1/N)

wy =2*pi*(0:(M-1))/M; %Create unshifted default omega axis
fy = 1/ps*(wy-pi*(1-mod(N,2)/N))/2/pi; %Shift zero to centre - for even case, pull back by pi, for odd case by pi(1-1/N)
[Fx,Fy] = meshgrid(fx,fy);


Del2_Psi_xy = func; %(-1*k*(I2-I1)/dz) %Epsilon added to denominator for divide by zero exception
Del2_Psi_uv = fftshift(fft2(Del2_Psi_xy));
Psi_uv = Del2_Psi_uv./(-4*pi^2*(Fx.^2+Fy.^2+epsilon));
Psi_xy = ifft2(ifftshift(Psi_uv),'symmetric');
end