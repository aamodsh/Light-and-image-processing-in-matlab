function E0 = initfield(N)
%function E0 = initfield(N)
%Returns desired field at input plane
%%
if mod(N,2)==1
    M = N+1; %For odd N; for even N this should be N
else 
    M = N;
end


%Creates N by N matrix filled with zeros
% global E0;
% E0 = zeros(N,N);

%%Point Spread
%E0(N/2,N/2)=1;


%Square hole phase object
E0 = ones(N,N);%E0(20:22,20:22)=exp(1i*pi);
size = 6 ; %Square size = 2*size  pixels
for index = (M/2-size):(M/2+size)
    for jindex = (M/2-size):(M/2+size)
        E0(index,jindex) = 0.25*exp(1i*pi);             %Gaussian phase exp(1i*(index^2+jindex^2)*2*pi);
     end
end


%90deg edges like litho mask
start = 42;
sizex = 1 ; %Square size = size  pixels
sizey = 6;
for index = start : (start+sizex)
    for jindex = (M/2-sizey):(M/2+sizey)
        E0(jindex,index) = exp(1i*pi/2);             %Gaussian phase exp(1i*(index^2+jindex^2)*2*pi);
     end
end


%90deg edges like litho mask
start = 54;
sizex = 1 ; %Square size = size  pixels
sizey = 6;
for index = start : (start+sizex)
    for jindex = (M/2-sizey):(M/2+sizey)
        E0(jindex,index) = exp(1i*pi/2);             %Gaussian phase exp(1i*(index^2+jindex^2)*2*pi);
     end
end





% %Single Slit
%   slitsize = 80; %in pixels
%   slitheight = 20;
  
%   %Random phase at each point
%   x=(M/2-slitheight/2):(M/2+slitheight/2);
%   y = (M/2-slitsize/2):(M/2+slitsize/2);
%  
%   for xind=1:length(x)
%       for yind = 1:length(y)
%           E0(xind,yind) = exp(1i*(rand(1))*2*pi);
%       end
%   end

% figure(10);imagesc(angle(E0));figure

% %Double Slit
% slitsize = N-2; %in pixels
% ypos = 50; %in pixels
% % E0((M/2-slitsize/2):(M/2+slitsize/2),(M/2-ypos)) = 1; %%Vertical
% % E0((M/2-slitsize/2):(M/2+slitsize/2),(M/2+ypos)) = 1;
% E0((M/2-ypos),(M/2-slitsize/2):(M/2+slitsize/2)) = 1; %%Horizontal
% E0((M/2+ypos),(M/2-slitsize/2):(M/2+slitsize/2)) = 1;


% % Sinusoidal
% period = 10 ; 
% x= linspace(-M/2,M/2,N);
% for i = 1:N
% E0(1:N,i) = sin(2*pi/period*x);
% end

% %Square periodic
% period = 50;
% for i = 1:N
%     if mod(floor(i./period),2)==1
%     E0(i,:)=1;
%     end
    

%figure(3);imagesc(abs(E0).^2);colormap(gray);
end
