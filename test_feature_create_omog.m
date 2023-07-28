%PICK OBJECT

%Lena-barbara
% E0 = sqrt(lena/255);
% Phi = barbara/255;
% test_feature = E0.*exp(1i*Phi);

%Shifted blocks
% test_feature = ones(568);Phi = ones(568);E0 = ones(568);
% Phi(200:300,200:300) = pi/6;
% E0(250:350,250:350) = 0.25;
% test_feature = E0.*exp(1i*Phi);


%Boundary layers
test_feature = ones(568);
absorber_angle = 0; %In degrees
absorber_alpha = 0.14; %alpha value is field value
bl_angle = 90; %across 6 pixel before gaussian spreading.

% test_feature(30:70,30:70) = absorber_alpha*exp(1i*absorber_angle*pi/180);
% test_feature(30:70,30:40) = 1*absorber_alpha*exp(1i*absorber_angle*pi/180)+1*exp(1i*bl_angle*pi/180); %add boundary layers
% test_feature(30:70,60:70) = 1*absorber_alpha*exp(1i*absorber_angle*pi/180)+1*exp(1i*bl_angle*pi/180);

%%Boundary layers
test_feature(254:304,254:304) = absorber_alpha*exp(1i*absorber_angle*pi/180);
test_feature(254:304,254:259) = 1*absorber_alpha*exp(1i*absorber_angle*pi/180)+1*exp(1i*bl_angle*pi/180); %add boundary layers
test_feature(254:304,298:303) = 1*absorber_alpha*exp(1i*absorber_angle*pi/180)+1*exp(1i*bl_angle*pi/180);


%%Cross amp phase
% test_feature(204:354,254:304) = absorber_alpha*exp(1i*absorber_angle*pi/180);
% test_feature(254:304,204:354) = absorber_alpha*exp(1i*bl_angle*pi/180); %add boundary layers
% test_feature(254:304,298:303) = 1*absorber_alpha*exp(1i*absorber_angle*pi/180)+1*exp(1i*bl_angle*pi/180);


%Phantom
% test_feature = 1.*exp(1i.*phantom(568));


%Phase vortex
% test_feature(254:304,254:304) = absorber_alpha*exp(1i*absorber_angle*pi/180);
% 
% x= 1:568;
% [X,Y] = meshgrid(x,x);
% temp_mat = unwrap(atan((Y-284)./(X-284)))/2/pi;
% figure;imagesc((temp_mat));colormap gray;
% 
% test_feature = test_feature.*exp(1i.*temp_mat);

%Square phase
% test_feature = ones(568);
% absorber_angle = 30; %In degrees
% absorber_alpha = 0.14; %alpha value is field value
% bl_angle = 90; %across 6 pixel before gaussian spreading.
% test_feature(254:304,254:304) = 1*exp(1i*absorber_angle*pi/180);
% test_feature(254:304,254:304) = 1*exp(1i*absorber_angle*pi/180);


% % %SMOOTHEN
smooth_sigma = 10; % is nominal
smooth_filt = fspecial('gaussian',[100 100],smooth_sigma); %Smoothen function, 30 pixels from lambda/NA
smooth_filt = smooth_filt/sum(smooth_filt(:));
Atemp = imfilter(abs(test_feature),smooth_filt,'replicate');
Phi_temp = imfilter(angle(test_feature),smooth_filt,'replicate');
test_feature = Atemp.*exp(1i*Phi_temp);


%SIMULATE DEFOCUS
figure; subplot(1,2,1);imagesc(abs(test_feature).^2);title('feature mag');axis off;colorbar;subplot(1,2,2);imagesc(angle(test_feature));title('feature phase');colorbar;axis off;




%%
%Initialize with/without noise
% noise_sig = 30/10000; %Experimental sigma = 30/10000
% I1=abs(test_feature_minus).^2;% + noise_sig*randn(568);
% I0=abs(test_feature).^2;%+ noise_sig*randn(568);
% I2=abs(test_feature_plus).^2;%+ noise_sig*randn(568);
% %
% I1=iso_direct_att_x_pol_reg(:,:,5);
% I0= iso_direct_att_x_pol_reg(:,:,6);
% I2=iso_direct_att_x_pol_reg(:,:,7);

% figure;imagesc(I2-I1);colorbar;title('Intensity difference');


% %First guess at recovered phase
% [phase_rec,phi_xy] = tie(I1,I0,I2,ps,lambda,2*dz,epsilon1,epsilonI);
% figure;imagesc(phase_rec);title('TIE phase recovered');colorbar;

