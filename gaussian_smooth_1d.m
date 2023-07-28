function outmat = gaussian_smooth_1d(test_feature,smooth_sigma)
%outmat = gaussian_smooth(test_feature,smooth_sigma)
%Smooths the amplitude and phase of complex field
smooth_filt = fspecial('gaussian',[101 101],smooth_sigma); %Smoothen function, 30 pixels from lambda/NA
smooth_filt = smooth_filt(51,:)/sum(smooth_filt(51,:));
Atemp = imfilter((test_feature),smooth_filt,'replicate');
outmat = Atemp;
end