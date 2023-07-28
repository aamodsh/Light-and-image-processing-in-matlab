stack = zeros(500,500,10);
stackf = zeros(500,500,10);

for q = 1:10

    ps = 3.6e-06;        %pixel size
    lambda = 632e-9;     %wavelength

    %% create gaussian
    sigma = 1.2;
    x = -249:250;
    [X,Y] = meshgrid(x,x);
    norm2D = (1/(2*pi))*exp(-0.5*(1/sigma^2).*(X.^2 + Y.^2));
    norm2D = norm2D - min(min(norm2D));
    norm2D = norm2D./(max(max(norm2D)));

    %% random phase object
    amp = zeros(500); 
    amp(101:400,101:400) = 1;
    phi = randn(500);
    phi = ifft2(fft2(phi).*fft2(norm2D));

    phi = phi./(max(max(abs(phi))))*pi;
    field = exp(1i.*phi); 

    %% propagate
    propf = field;
    for i = 1:10
        propf = propagate(propf,lambda,10e-4,ps,500);
        stack(:,:,i) = stack(:,:,i)+abs(propf).^2;
    end
end

stack = stack./q;

%% plots

for p=1:10
    I = stack(:,:,p);
    IF = abs(fftshift(fft2(I)));
    IF(251,251) = 0;
    stackf(:,:,p) = IF;
    figure(p); imagesc(IF)
    pause(0.3)
end
