%MULTIPLE IMAGE TIE

N=800; %Input from splat at wafer plane
dz = 100;
epsilon = 0.8*10^-8;

stack_size = 5; %No of input images



Ivec = zeros(N,N,stack_size);
xvector = linspace(-(stack_size-1)/2,(stack_size-1)/2,stack_size);
dzvec = dz*xvector;

%%
% E0 = initfield(N);
% %generate images 
% for i = 1:stack_size
% Evec(:,:,i) = fresnel_prop(E0,N,30,193,dzvec(i));
% end
% 
% Evec(:,:,(stack_size+1)/2)=E0;


% % %Read Images
% Ivec(:,:,1)=double(I_1);
% Ivec(:,:,2)=double(I_1);
% Ivec(:,:,3)=double(I_2);




% % % %Read Images
Ivec(:,:,1)=M_200(100:899,300:1099);
Ivec(:,:,2)=M_100(100:899,300:1099);
Ivec(:,:,3)=M0;
Ivec(:,:,4)=M100(100:899,300:1099);
Ivec(:,:,5)=M200(100:899,300:1099);


% % %Read Images
% Ivec(:,:,1)=A_8;
% Ivec(:,:,2)=A_6;
% Ivec(:,:,3)=A_4;
% Ivec(:,:,4)=A_2;
% Ivec(:,:,5)=A0;
% Ivec(:,:,6)=A2;
% Ivec(:,:,7)=A4;
% Ivec(:,:,8)=A6;
% Ivec(:,:,9)=A8;





%Run TIE from intensities
[Phi,Psi,gradx,grady] = tie_stack((Ivec),30,193,dz,epsilon);
figure(2);imagesc(Phi);colorbar;colormap gray;
% rectangle('Position',[43 43 12 12],'Linewidth',2,'EdgeColor','black');

% %Run TIE from fields
% [Phi,Psi,gradx,grady] = tie_stack(abs(Evec).^2,30,193,dz,epsilon);
% figure;imagesc((Phi));colormap gray;


%%
%%Run tie for different values of epsilon
% exponent=linspace(2,10,4);
% %Parametrize with epsilon
% figure(6);
% 
% for i = 1:4; 
% epsilon=10^exponent(i);%(10+(2*i/12));
% [Phi,Psi,gradx,grady] = tie_stack(Evec.^2,30e-9,193e-9,dz,epsilon);
% subplot(2,2,i);
% imagesc(real(Phi));colorbar; caxis([-pi pi]);title(num2str(exponent(i)));
% 
% end