%
%SINGLE dz TIE

dz=100;
N=97*6 ;


%Create mask and progpagated fields
E0 = initfield(N); %Create mask replica in matlab
E1 = fresnel_prop(E0,N,30,193,-1*dz);
E2 = fresnel_prop(E0,N,30,193,dz);

% %Turn to intensity
% I0 = abs(E0).^2;
% I1 = abs(E1).^2;
% I2 = abs(E2).^2;
 
%Add poisson noise to intensities
% n=0.1;
% I0 = imnoise(n*(abs(E0).^2),'poisson');
% I1 = imnoise(n*(abs(E1).^2),'poisson');
% I2 = imnoise(n*(abs(E2).^2),'poisson');

%Moustache man measured object
% I0 = M0;
% I1 = M_100(100:899,300:1099);
% I2 = M100(100:899,300:1099);





% figure;imagesc(unwrap(angle(E0)));colormap gray;title('object');

%%
% %One step tie
% figure;
% epsilon = 1e-25;epsilonI = 0;
% [Phi,Psi,gradx,grady] = tie((I1),(I0),(I2),30,193,2*dz,epsilon,epsilonI);
% im = imagesc((unwrap(Phi)));colorbar; colormap gray; %caxis([-pi pi]);
% title(num2str(epsilon));
% error = sum(sum((unwrap(Phi)-unwrap(angle(E0)).^2)))/N/N;

%%
%TIE with noise parametrization
% figure;
% levels = 8;
% n_photon = linspace(1,8,levels);
% error = zeros(levels,1);
% 
% for nindex = 1:levels
%     %Add poisson noise to intensities
% 
%     I0 = imnoise(n_photon(nindex)*(abs(E0).^2),'poisson');
%     I1 = imnoise(n_photon(nindex)*(abs(E1).^2),'poisson');
%     I2 = imnoise(n_photon(nindex)*(abs(E2).^2),'poisson');
% 
%     epsilon = 1e-10;epsilonI = 0;
%     [Phi,Psi,gradx,grady] = tie((I1),(I0),(I2),30,193,2*dz,epsilon,epsilonI);
% 
%     subplot(levels/4,4,nindex); 
%     im = imagesc((unwrap(Phi)));colorbar; colormap gray; %caxis([-pi pi]);
%     title(num2str(n_photon(nindex)));
%     error(nindex) = sum(sum((unwrap(Phi)-unwrap(angle(E0)).^2)))/N/N;
% end
% 
% figure;plot(n_photon,error);

%%
%TIE with epsilon parametrization
% Vary exponent value over loop no of cycles
% loops=32;
% epsilon_vec=10.^linspace(-12,-3,loops);
% error = zeros(loops,1);
% 
% % Parametrize with epsilon
% figure(6);
% %Run tie for different values of epsilon
% for i = 1:loops; 
% epsilon=epsilon_vec(i);%(10+(2*i/12));
% [Phi,Psi,gradx,grady] = tie(abs(I1),abs(I0),abs(I2),30,193,2*dz,epsilon,0);
% subplot(loops/4,4,i);
% imagesc(unwrap(Phi));colormap gray;colorbar; %caxis([-pi pi]);
% % title(num2str(epsilon_vec(i)));
% title(strcat(num2str(epsilon_vec(i)),'; photons=',num2str(n)));
% 
% error(i) = sum(sum((unwrap(Phi)-unwrap(angle(E0))).^2))/N/N;
% 
% end
% 
% 
% %plot mean square the error with expected phase
% figure(5);hold on;
% plot(log10(epsilon_vec),error);

%%
% TIE variation with epsiLon with noise as parameter 
nlevels = 30; %No of parameter steps
loops =32; %No of epsilon steps to loop over for each plot (multiple of 4 ideally for subplot)
n_photon = linspace(0.1,1.5,nlevels); %Parameter range for noise
Leg = strread(num2str(n_photon),'%s');%Legend
error = zeros(loops,nlevels);
epsilon_vec=10.^(linspace(-12,-3,loops)); %Parameter range for regularization

% Vary noise parameter
for nindex = 1:nlevels
%     %Add poisson noise to intensities
% 
    I0 = imnoise(n_photon(nindex)*abs(E0).^2,'poisson');
    I1 = imnoise(n_photon(nindex)*abs(E1).^2,'poisson');
    I2 = imnoise(n_photon(nindex)*abs(E2).^2,'poisson');
    
     
    % Parametrize with epsilon
%     figure(nindex);

    %Run tie for different values of epsilon
    for i = 1:loops; 
        epsilon=epsilon_vec(i);%(10+(2*i/12));
        [Phi,Psi,gradx,grady] = tie(abs(I1),abs(I0),abs(I2),30,193,2*dz,epsilon,0);
%         subplot(loops/4,4,i);
%         imagesc(unwrap(Phi));colormap gray;colorbar; %caxis([-pi pi]);
%         title(strcat(num2str(epsilon_vec(i)),'; photons=',num2str(n_photon(nindex))));
        error(i,nindex) = sum(sum((unwrap(Phi)-unwrap(angle(E0))).^2))/N/N; %Avg MS error per pixel

    end
end
    
figure(5); plot(log10(epsilon_vec),error);
legend(Leg);


%TIE with epsilon parametrization and no noise

%Intensity without noise
I0 = abs(E0).^2;
I1 = abs(E1).^2;
I2 = abs(E2).^2;


% Vary exponent value over loop no of cycles
error0 = zeros(loops,1);

% Parametrize with epsilon
% figure;
%Run tie for different values of epsilon
for i = 1:loops; 
epsilon=epsilon_vec(i);%(10+(2*i/12));
[Phi,Psi,gradx,grady] = tie(abs(I1),abs(I0),abs(I2),30,193,2*dz,epsilon,0);
% subplot(loops/4,4,i);
% imagesc(unwrap(Phi));colormap gray;colorbar; %caxis([-pi pi]);
% title(num2str(epsilon_vec(i)));
% title(strcat(num2str(epsilon_vec(i)),'; photons=',num2str(n)));

error0(i) = sum(sum((unwrap(Phi)-unwrap(angle(E0))).^2))/N/N;

end


%plot mean square the error with expected phase
figure(5);hold on;
plot(log10(epsilon_vec),error0,'--k');


%%
%TIE with dz PARAMETRIZATION
%Error versus dz
% dzvec = 1:40:310;
% error = zeros(length(dzvec));
% epsilon = 10^-10;
% N=97*6;
% 
% E0 = initfield(N);I0 = abs(E0).^2;
% 
% for index = 1:length(dzvec) 
%     dz = dzvec(index);
%     
%     E1 = fresnel_prop(E0,N,30,193,-1*dz);
%     E2 = fresnel_prop(E0,N,30,193,dz);
%     I1=abs(E1).^2;
%     I2=abs(E2).^2;
%  
%     [Phi,Psi,gradx,grady] = tie((I1),(I0),(I2),30,193,2*dz,epsilon,0);
%     subplot(length(dzvec)/4,4,index);imagesc(Phi);title(strcat('z=',num2str(dz)))
%     error(index) = sum(sum((unwrap(Phi-unwrap(angle(E0))).^2)))/N/N;
% end
% 
% figure;plot(dzvec,error);




