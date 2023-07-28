function Im0 = phase_object(N)
M=200;
phase = phantom(M);
phase = padarray(phase,[(N-M)/2 (N-M)/2]); 

%Gaussian to smooth
G = fspecial('gaussian',[N N],N/400);
% imagesc(G);
phase2 = imfilter(phase,G,'same');

Im0 = ones(N).*exp(1i*2*pi*phase2);
% figure;imagesc(angle(phaseobj));colormap gray;
