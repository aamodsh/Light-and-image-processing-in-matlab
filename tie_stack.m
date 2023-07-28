function [Phi_xy,Psi_xy,Grad_Psi_x,Grad_Psi_y,grad2x,grad2y] = tie_stack(Ivec,ps,lambda,dz,epsilon)
%Accepts an odd number of 2d images as 3D matrix Ivec, the stack of images
%generated from a propagator,and the z distance between planes :      dz
%The images are centred about the plane being solved for.
%Generates the phase profile from the inverse laplacian, and back
%propogates. The guessed imaged is compared with original image:    Image 
%as a mean square error. 
%Pixel Size :                                                       ps(1um)
%Wavelength :                                                       lambda(0.628)
%Ouput      :                                                       Reconstructed Phase 

%I0 = I1; %Central plane; can be passed as an argument, here I0~I1
stack_size = size(Ivec,3);
I0 = Ivec(:,:,(stack_size+1)/2); %central plane
N = size(Ivec,1);
k = 2*pi/lambda;
xsize =  ps*(N-1); ysize = ps*(N-1); %square grid with N*N points
%epsilon = 01e1;

%Create space
x = linspace(-xsize/2,xsize/2,N);   %N point sampling over xsize
y = linspace(-ysize/2,ysize/2,N);

%Inverse space
% if mod(N,2)==1
%     M = N+1; %For odd N; 
% else 
    M = N;
% end

%Proper way of creating frequency axis
wx =2*pi*(0:(N-1))/N; %Create unshifted default omega axis
%fx =1/ps*unwrap(fftshift(wx)-2*pi)/2/pi; 
fx = 1/ps*(wx-pi*(1-mod(N,2)/N))/2/pi; %Shift zero to centre - for even case, pull back by pi, for odd case by pi(1-1/N)
[Fx,Fy] = meshgrid(fx,fx);

%%
%Other way
% fx = linspace(-(N-1)/2*(1/xsize),(N-1)/2*(1/xsize),N);%Fx = repmat(fx,N,1);
% fy = linspace(-(N-1)/2*(1/ysize),(N-1)/2*(1/ysize),N);%Fy = repmat(fy',1,N);
%[Fx,Fy] = meshgrid(fx,fy);
%EXPLANATION - original fft goes from 1-N in N points.
%Hence (N-1)/2 points on either side of zero = 2*(N-1)/2 + 1 = N points total
% Fx(Fx==0) = 1; Fy(Fy==0)=1;                 %Important to prevent Inf/NaN from appearing

% [val,index_centre] = min(abs(fx)); %Finds the index of centre element (useful in plotting)


%%
%Calculate the derivative using all the stacks 

%Calculate finite different coefficients
xvector = linspace(-(stack_size-1)/2,(stack_size-1)/2,stack_size);
c = TT(xvector,1)/dz; %Centred about 0, c is the taylor coefficients (predivide by dz)

% 
% tic;
% for i = 1:length(c) 
% Cmat(:,:,i) = repmat(c(i),N,N);
% end
% toc;


%Apply coefficients
% taylor_mat = Cmat.*Ivec;
% dIdz = sum(taylor_mat,3); %Add all matrices along third dimension

%Sam's fast fix
dIdz = zeros(size(Ivec,1),size(Ivec,2));
for i=1:length(c)
    dIdz = dIdz + c(i)*squeeze(Ivec(:,:,i));
end

%%
%Solve the Lapacian
%RHS of laplacion del2(I.Del(phi)) = Del2_Psi
% Del2_Psi_xy = (-1*k*(I2-I1)/dz)%/(I0+1e-6); %Epsilon added to denominator for divide by zero exception
% Del2_Psi_uv = fftshift(fft2(Del2_Psi_xy));
% Psi_uv = Del2_Psi_uv./(-4*pi^2*(Fx.^2+Fy.^2+1e-6));
% Psi_xy = ifft2(ifftshift(Phi_uv));
Del2_Psi_xy = (-1*k*dIdz); %Obtain the derivative from stack
Psi_xy = poisson_solve(Del2_Psi_xy,ps,N,epsilon);
[Grad_Psi_x, Grad_Psi_y] = gradient(Psi_xy/ps); %Take the x and y gradients
Grad_Psi_x = Grad_Psi_x./I0;Grad_Psi_y = Grad_Psi_y./I0; %Divide by intensity
[grad2x,dummy1] = gradient(Grad_Psi_x/ps);[dummy2,grad2y]=gradient(Grad_Psi_y/ps);
Del2_Phi_xy = grad2x +grad2y;
Phi_xy = poisson_solve_symm(Del2_Phi_xy,ps,N,epsilon);

%Renormalize (till the scaling factor is figured out)
%   maxval=max(max((real(Phi_xy))));
%   minval=min(min((real(Phi_xy))));
%   Phi_xy = (Phi_xy-minval)/(maxval-minval)*pi;
%   Phi_xy = angle(exp(1i.*Phi_xy)); 


%Substract dc term and multiply by minus (emperical)
Phi_xy = -1*(Phi_xy - sum(Phi_xy(:,1))/length(Phi_xy(:,1)));
%%Wrap angle
%Phi_xy = angle(exp(1i*Phi_xy));



end