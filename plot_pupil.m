function [kx_vec ky_vec puprange_ext NA_min NA_max] = plot_pupil_kx(NA,N,lambda,ps)
%plot_pupil(NA,N,lambda,ps)
%Returns 
%kx_mat_puprange : pupil values within the NA  
%puprange_ext : range of pixel# within the NA

full_NA_pixels = kaxes(lambda,ps,N); %depends on ps,lambda,N
NA_range = (full_NA_pixels-1)/2;  %number of samples on either side of zero for NA =1
                                  %obtained using kaxes.m (with NA =1)

mid = (N+2)/2;
upper_index = round(mid + round(NA_range*NA));
lower_index = round(mid - round(NA_range*NA));
pupil = zeros(1,N);
pupil(lower_index:upper_index) = 1;

%%--plot the pupil a little outside the NA as well
puprange = find(pupil~=0); %spans the whole pupil for NA_in = 1;
pupextra = round(0*length(puprange)); %percentage of hyperNA to display - 5% by default
puprange_ext = (puprange(1)-pupextra):((puprange(length(puprange))+pupextra));


%real pupil coordinates
kpupil = (1:length(puprange))/ps/N;
[Xtemp,Ytemp] = meshgrid((1:N)-(N/2+1),(1:N)-(N/2+1));
kx_mat = Xtemp/N/ps*2*pi;
NA_min = kx_mat(1,puprange(1));
NA_max = kx_mat(1,puprange(length(puprange)));

kx_mat_puprange = kx_mat(1,puprange);
kx_vec = kx_mat(1,:);

NA_min = kx_mat(1,puprange(1));
NA_max = kx_mat(1,puprange(length(puprange)));


end

                                  
                                  