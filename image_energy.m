function im_E = image_energy(images)
N_frame = size(images,3);
im_E = zeros(N_frame,1);
for k = 1:N_frame
    im_line = images(:,:,k);
    im_line = im_line(:);
    im_E(k) = sqrt(sum(im_line.^2));
end
end