function [Im2reg,xshift,yshift] = register_ims(Im1,Im2)
%function [Im2reg] = register_ims(Im1,Im2)
%Register image 2 to image 1; output registered image/shift values
%% NOTE: some problems with Im2reg, the min values are not preserved as for Im2
[output Greg] = dftregistration(fft2(Im1),fft2(Im2),30); %Third argument is sub-pixel discretization
yshift = output(3)
xshift = output(4)
% Im2reg = subpixel_shift(Im2,yshift,xshift);
Im2reg = (ifft2(Greg,'symmetric'));
end