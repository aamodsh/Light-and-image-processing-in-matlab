function [outint outphase Efftpupil pupil] = scatter_NA_filter(nfmag,nfphase,c,sigma,NA_in,lambda,ps)
%scatter_NA_filter(nfmag,nfphase,c)
%%Aerial imaging with low pass filtering and partial coherence,
%%(c)Computational imaging lab 2015 (aamod@berkeley.edu)
%Inputs are mask nearfields and coherence (c = 1 is coherent, c = 0 is partially coherent with illumination sigma as entered)
%outputs filtered intensity, phase, and full pupil
%% Mask near field
%for NA=1
N = length(nfmag);
full_NA_pixels = kaxes(lambda,ps,N); %depends on ps,lambda,N
NA_range = (full_NA_pixels-1)/2;  %number of samples on either side of zero for NA =1
                                  %obtained using kaxes.m (with NA =1)
%Imaging system NA
NA = NA_in; %input from function call
Efield = nfmag.*exp(1i.*nfphase);

%---Parital coherence filter---
% sigma = 0.3;
klim_illu = round(sigma*NA_range*NA); %kspace edge pixel of illumination
if c==0
    nsrc = 7;
    shift = -klim_illu:klim_illu; %Crude partial coherence (pixels(2+2+1), 2 pixels = 20*0.3375*0.3)
else
    shift = [0]; %Coherent
end


%---Implement imaging---%
clear Efft Efftpupil Efftpupil3D Efield_wafer

for i = 1:length(shift)
    
Efft = fftshift(fft(Efield));
Efft_shift = circshift(Efft,[0 shift(i)]);
%* Efft_shift = subpixel_shift_1d(Efft,shift(i));

mid = (length(Efft)+2)/2;
upper_index = round(mid + round(NA_range*NA));
lower_index = round(mid - round(NA_range*NA));
pupil = zeros(1,length(Efft));
pupil(lower_index:upper_index) = 1;

%Smooth pupil
% Efftpupil = Efft.*gaussian_smooth_1d(pupil,1);
%Coherent pupil
% Efftpupil = Efft.*pupil;

%Return pupil for each source point (3D pupil)
Efftpupil3D(i,:) = Efft_shift.*pupil;
%Or return coherent pupil
Efftpupil = Efft;


% figure;plot(abs(Efftpupil(i,:))); title(['Pupil for source' num2str(i)]);
Efield_wafer(i,:) = ((ifft(ifftshift(Efftpupil3D(i,:)))));

end

%Recover wafer
%Output field magnitude
outint = sum(abs(Efield_wafer(:,:)).^2,1);

% %Output field phase
midshift = (length(shift)+1)/2;
if c == 0
 outphase = angle(Efield_wafer(midshift,:));
    for i  = 1:(length(shift)-1)/2
        outphase = outphase + unwrap(angle(Efield_wafer(midshift+i,:)))+...
            unwrap(angle(Efield_wafer(midshift-i,:)));
    end
outphase = outphase/length(shift);
else 
    outphase = (angle(Efield_wafer(1,:)));
end
