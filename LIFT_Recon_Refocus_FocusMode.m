% ***** 
% The main script to do focus mode reconstruction and refocus for depth-from-focus
% *****
close all; clear all; clc;
topDIR = '.\LIF\DeSense_Ver\';
data_DIR = [topDIR 'SpotPattern\'];  % CalibData\Depth\
reconFile = 'SpotPattern';  
lambda = 1e-1;    % the regularization parameter for FISTA

% 1: Load the calibration data results
Calib_Res = load([topDIR 'Calib_Res.mat']);
Calib_Res = Calib_Res.Calib_Res;
image_t = load([data_DIR reconFile '.mat'], 'image_t');
image_t = image_t.image_t;
image_t = norm1(image_t); 

%% *****
% shearing and reconstruction
% ******
options.INVERT = false;    options.CROP = true;   options.DECONV = true;
options.Refocus = true; 
options.USE_TV = false;  % whether to use TV or l1, l1 is default
% refocusing by changing the lenslet center: shearing light field
s_array = linspace(1,-2,30); % shearing factor
u = (1:7)-4;  % angular resolution 
Nstep = length(s_array);

cntx_depth = Calib_Res.cntx_depth;
sub_img_cnt_new = round(cntx_depth(:,4)).'; % cntx_depth index differet nominal image plane 
start_t = tic;
for K = 1:Nstep
    sub_img_cnt = round(sub_img_cnt_new + u.*s_array(K));
    options.sub_img_cnt = sub_img_cnt;
    cnt_k(:,K) = sub_img_cnt;
    % 3: divide the slice into N_angle lenslet sub-regions
    im_crop(:,:,K) = fx_LIFT_Recon2D(Calib_Res, image_t, lambda, options);
end
disp('shearing and reconstruction time:');
run_t = toc(start_t)

figure('position', [200, 200, 1600, 600]) 
montage(norm1(im_crop),'DisplayRange', [0 1]); colormap('hot');

%% ***** 
%  perform depth from focus
%  *****
% 1: align the image focus stack
% use the first image as the reference image
[optimizer, metric] = imregconfig('multimodal');
tform_Cell = cell(Nstep,1);
for K = 2:size(im_crop,3)
    im_ref = squeeze(im_crop(:,:,K-1));
    % tform_Cell{K} = imregcorr(im_crop(:,:,K),im_ref); 
    tform_Cell{K} = imregtform(im_crop(:,:,K), im_ref, 'similarity', optimizer, metric); 
end

Rfixed = imref2d(size(im_ref));
for K = 2:size(im_crop,3)
    if(K==2)
        tform = tform_Cell{2};
    else
        tform = affine2d(tform.T * tform_Cell{K}.T);
    end
    im_crop(:,:,K) = imwarp(im_crop(:,:,K),tform,'OutputView',Rfixed); 
end

% 2 focal stack filtering: vbm3d 
im_crop = norm1(im_crop); % normalize across the whole focal stack is fine
[PSNR, im_crop] = VBM3D(im_crop, 50);
im_crop(isnan(im_crop)) = 0.0;

% refocusd focal stack visualization
figure('position', [200, 200, 1600, 600]) 
montage(padarray(norm1(im_crop),[1,2],1),'DisplayRange', [0 1], 'size',[6,5]); colormap('hot');

% 3: DfF calculation using focal measure
im_crop_stack = im_crop;
Method = 'SML';  % Focus measure: selecting from {'IME', 'SPARC', 'RSML', 'SML'}
[depth_fit, depth_map, Img_ALLFOCUS] = depth_from_focus(im_crop_stack, Method, 0.05, false, false);

figure; subplot(1,2,1); imagesc(depth_fit); colormap('hot');
subplot(1,2,2); surf(depth_fit.','LineStyle','None'); colormap('hot'); title('Fit depth surface');

depth_fit = medfilt2(depth_fit,[5,5]);

% **********
% convert idx of shearing to physical depth: not needed for
% relative depth extraction
% **********
a_opt = Calib_Res.b_opt;
pixelSize = Calib_Res.pixelSize;
delta_depth_est_Func = @(dx) a_opt(1).* a_opt(3)./(dx*pixelSize-a_opt(3));
% the baseline/reference plane induced RegSize
base_dist = mean(sub_img_cnt_new(2:end) - sub_img_cnt_new(1:end-1)); 
dist_refocus = base_dist + s_array;
depth_estimate = delta_depth_est_Func(dist_refocus);
depth_fit(depth_fit<1) =1;
depth_phy = depth_estimate(depth_fit);
depth_phy = depth_phy - min(depth_phy(:));

figure; surf(depth_phy,'LineStyle','None'); colormap('hot'); 
axis square; title('3D Depth map'); axis([0,129,0,129]); view(3);
set(gca,'FontSize',18)