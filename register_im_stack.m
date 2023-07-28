function [stack_reg xshift yshift] = register_im_stack(im_stack)
%function stack_reg = register_im_stack(im_stack)
%registers image stack using neighbouring images in sequence
    stack_reg(:,:,1) = im_stack(:,:,1);
    for i = 2:size(im_stack,3)
       [stack_reg(:,:,i) xshift(i) yshift(i)] = register_ims(stack_reg(:,:,i-1),im_stack(:,:,i)); 
    end

