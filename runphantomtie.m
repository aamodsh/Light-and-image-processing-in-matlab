dz=100;
N=512;
Im1 = fresnel_prop(initfield(N),N,30,193,-1*dz);
Im2 = fresnel_prop(initfield(N),N,30,193,dz);
E0  = (phase_object(N));
figure;imagesc(unwrap(angle(E0)));colormap gray;title('object');

%%
% % %TIE with epsilon parametrization
% % %Vary exponent value over loop no of cycles
% loops=16;
% epsilon_vec=10.^(linspace(-15,1,loops));
% error = zeros(loops,1);
% %Parametrize with epsilon
% figure(6);
% 
% %Run tie for different values of epsilon
% for i = 1:loops; 
% epsilon=epsilon_vec(i);%(10+(2*i/12));
% [Phi,Psi,gradx,grady] = tie(abs(Im1).^2,abs(E0).^2,abs(Im2).^2,30,193,2*dz,epsilon);
% subplot(loops/4,4,i);
% %imagesc(Phi-angle(E0));colorbar; caxis([-pi pi]);title(num2str(exponent(i)))
% imagesc(unwrap(Phi));colormap gray;colorbar; %caxis([-pi pi]);
% title(num2str(epsilon_vec(i)));
% %
% error(i) = sum(sum((unwrap(Phi)-unwrap(angle(E0))).^2))/N/N;
% 
% end
% 
% 
% %plot mean square the error with expected phase
% figure;plot(log10(epsilon_vec),error);

%%
%One step tie
figure;
epsilon = 10^-10;
[Phi,Psi,gradx,grady] = tie(abs(Im1).^2,abs(E0).^2,abs(Im2).^2,30,193,2*dz,epsilon);
imagesc(unwrap(Phi));colorbar; colormap gray; %caxis([-pi pi]);
title(num2str(epsilon));


