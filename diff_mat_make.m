%naya diffuser banein
nlambda = 50;
for j = 1:nlambda
ps = 3.6e-6; %andor camera pixel size
lambda = linspace(1,1000,nlambda)*1e-9; %wavelength sweep

%%----create axis
% N = 900;
% 
% figure;imagesc(X)

%%--MAKE DIFFUSER--%
% roughness_ind = 1/4; %high roughness (amp/fwhm)
% sigma = 11; %standard deviation
% scale = sigma*roughness_ind;
% gaussian2D = scale*1/(2*pi*sqrt(sigma))*exp(-1*(X.^2 + Y.^2)./(2*(sigma).^2));
% % imagesc(X(1,:), Y(:,1), gaussian2D);
% % mesh(gaussian2D);axis([0 500 0 500 0 1]);
% 
% randn = rand(N,N);
% % figure;imagesc(randn);
% 
% %%convolve 
% diff_surf = ifft2((fft2(randn).*fft2(gaussian2D)));
% % figure;mesh(real(diff_surf));colormap(gray);
% 
% %%--diffuser_make 
%  prop_zero = exp(1i.*diff_surf); %in radian, consider normalizing by lambda
%  %%-------------------
 
 %% ---CREATE PROPAGATOR---%%
 nz = 33;
 zmin = -5e-4;
 zmax = 5e-4;
 
%defocus steps...
 zstack= linspace(zmin,zmax,nz); 
 clear propf  stack;
 for i = 1:nz
        propf = fresnel_prop(prop_zero, lambda(j),ps ,zstack(i) );
        stack(:,:,i) = abs(propf).^2;
%         
%         figure(14);
%         contour((stack(:,:,i)),25,'linewidth',1);
%         
%         figure(119+j);
        contour(sqrt(stack(:,:,i)),25,'linewidth',1);
        title(['lambda =' num2str(lambda(j)) 'z ='  num2str(zstack(i))]);
        %         f(:,:,j) = getframe;
%      pause(0.1);
%         lambda(j)
%       waitforbuttonpress;
end    
end 

%% save image
saveMovie(f,'strong_speckle_lambda_1_1000.gif',10,'k');
