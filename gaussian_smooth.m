function outmat = gaussian_smooth(test_feature,smooth_sigma)
%outmat = gaussian_smooth(test_feature,smooth_sigma)
%Smooths the amplitude and phase of complex field
smooth_filt = fspecial('gaussian',[100 100],smooth_sigma); %Smoothen function, 30 pixels from lambda/NA
smooth_filt = smooth_filt/sum(smooth_filt(:));
Atemp = imfilter(abs(test_feature),smooth_filt,'replicate');
Phi_temp = imfilter(angle(test_feature),smooth_filt,'replicate');
outmat = Atemp.*exp(1i*Phi_temp);
end