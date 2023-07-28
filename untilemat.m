function outmat = untilemat(inmat)
outmat = inmat(1:size(inmat,1)/2,1:size(inmat,2)/2);
end