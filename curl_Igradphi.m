% Given matrix I0 and phase phi, computes curl of IgradPhi = gradI curl
% gradPhi
function curlmat = curl_Igradphi(I0,Phi,ps)
%function curlmat = curl_Igradphi(I0,Phi,ps)

% phi = angle(test_feature);

[gradIx,gradIy] = gradient(I0/ps);
[gradphix,gradphiy] = gradient(Phi/ps);
curlmat = gradIx.*gradphiy-gradIy.*gradphix;

%Plot
% figure;subplot(2,2,1);imagesc(gradIx);colorbar;
% subplot(2,2,2);imagesc(gradIy);subplot(2,2,3);imagesc(gradphix);subplot(2,2,4);imagesc(gradphiy);colorbar;
% figure;imagesc(gradIx.*gradphiy);colorbar;
% figure;imagesc(gradIy.*gradphix);colorbar;
% figure;imagesc(curlmat);title('gradI x gradPhi');colorbar;axis image;axis off;
