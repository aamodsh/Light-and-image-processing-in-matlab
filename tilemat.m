function outputmat = tilemat(inputmat)

for i = 1:size(inputmat,3);
inmat = inputmat(:,:,i);
quad2 = fliplr(inmat);
quad3 = flipdim(quad2,1);
quad4 = flipud(inmat);
outmat = zeros(2*size(inmat,1),2*size(inmat,2));
outmat(1:size(outmat,1)/2,1:size(outmat,2)/2)=inmat;
outmat((size(outmat,1)/2+1):size(outmat,1),(size(outmat,2)/2+1):size(outmat,2))=quad3;
outmat((size(outmat,1)/2+1):size(outmat,1),1:size(outmat,2)/2)=quad4;
outmat(1:size(outmat,1)/2,(size(outmat,2)/2+1):size(outmat,2)) = quad2;
outputmat(:,:,i)=outmat;
end