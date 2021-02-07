% ************* This script perform LIFT reconstruciton in streaking mode/time domain *******************
% ****** Reconstruct the volume and perform refocusing for depth extraction ********
% ******  ***********************
close all; clear all; clc;

top_DIR = '.\LIF\DeSense_Ver\'; 
DIR = top_DIR; calib_DIR = top_DIR;
filename = 'HelicalFiber';
lambda_reg = 5e-1;   % the regularization parameter

% ******
% 1: Load the calibration data results figure
% ******
Calib_Res = load([calib_DIR,'Calib_Res.mat']);
Calib_Res = Calib_Res.Calib_Res;
streak_tform = Calib_Res.streak_tform;

% ******
% 2: Load image data and correct the streak distortion
% ******
image_t = load([DIR filename '\' filename '.mat']);
image_t = image_t.image_t;
figure(1); subplot(1,2,1); imagesc(image_t(1:end,:)); colormap('hot'); title('averaged image');

image_t = imwarp(image_t, streak_tform, 'OutputView', imref2d(size(image_t)));
subplot(1,2,2); imagesc(norm1(image_t)); title('rectified streak image');
drawnow();

% crop decimate the signal to reduce reconstruction time: depth from focus
% need texture/contrast: avoid blank region
image_t = image_t(100:3:720,:);

%% ***********
%  3: shearing and reconstruct
%  ***********
% setup solver
options = []; options.INVERT = false;  options.CROP = true; 
options.Deconv = true;  options.Normalize = false;
options.USE_TV = false;  % flag to indicate whether use TV as regularizer, default is l1
options.Refocus = true; % flag to indicate whether perform refocusing
% refocusing part: postive s --> closer, negative s --> further
s_array = linspace(4,-4,31);  % the shearing factor for refocusing
u = (1:7)-4; % the lenslet/angular resolution indexing
Nstep = length(s_array);
cntx_depth = Calib_Res.cntx_depth;
sub_img_cnt_new = round(cntx_depth(:,3)).';

start_t = tic;
for K = 1:Nstep
    % 1: shearing the light field: equivalent to shift the sub-image center
    sub_img_cnt = round(sub_img_cnt_new + u.*s_array(K));
    options.sub_img_cnt = sub_img_cnt;
    cnt_k(:,K) = sub_img_cnt;
    
    % 2: reconstruct
    im_crop(:,:,:,K) = fx_LIFT_ReconVOL(Calib_Res, image_t, lambda_reg, options);
end
run_t = toc(start_t)

%% ********
%  5: DFF: depth from focus method
% *********
im_frame_focus = zeros(size(im_crop,1), size(im_crop,1), size(im_crop,3));
idx_depth = zeros(size(im_crop,3),1);
Method = 'SML';  % Focus measure: selecting from {'IME', 'SPARC', 'RSML', 'SML'}
for K = 1:size(im_crop,3)  % the number of frames
    % normalize across the whole focal stack
    im_f_stack = norm1(squeeze(im_crop(:,:,K,:))); 
    
    [~,im_f_stack] = VBM3D(im_f_stack,10);
    im_f_stack(isnan(im_f_stack)) = 0.0;

    % For spatially-sparse scenes, image energy (L2 norm) can be an
    % effective focal meassure.
    im_E = image_energy(im_f_stack);
    im_E = medfilt1(im_E,13);
    [~, idx_depth(K)] = max(im_E);
    % All-in-focus by extracting the image parts with maximum focus measure
    im_frame_focus(:,:,K) = im_crop(:,:,K,idx_depth(K));    
end

%% convert idx of shearing to physical relative depth: optional step
a_opt = Calib_Res.b_opt;
pixelSize = Calib_Res.pixelSize;
delta_depth_est_Func = @(dx) a_opt(1).* a_opt(3)./(dx*pixelSize-a_opt(3));
% the baseline/reference plane induced RegSize
base_dist = mean(sub_img_cnt_new(2:end) - sub_img_cnt_new(1:end-1)); 
dist_refocus = base_dist + s_array;
depth_estimate = delta_depth_est_Func(dist_refocus);
depth_vol = depth_estimate(idx_depth);

%% all-in focus
im_all_in_focus_DC = sum(im_frame_focus,3);
figure; 
subplot(2,2,1); imagesc(sum(squeeze(im_crop(:,:,:,1)),3)); axis square; 
title('Refocus at 1'); colormap('hot'); set(gca,'XTick',[], 'YTick', [])
subplot(2,2,2); imagesc(sum(squeeze(im_crop(:,:,:,10)),3)); axis square;
title('Refocus at 10'); colormap('hot'); set(gca,'XTick',[], 'YTick', [])
subplot(2,2,3); imagesc(sum(squeeze(im_crop(:,:,:,15)),3)); axis square;
title('Refocus at 15'); colormap('hot'); set(gca,'XTick',[], 'YTick', [])
subplot(2,2,4); imagesc(sum(squeeze(im_crop(:,:,:,30)),3)); axis square; 
title('Refocus at 30'); colormap('hot'); set(gca,'XTick',[], 'YTick', []);
figure; imagesc(im_all_in_focus_DC); axis square; 
title('All in focus'); colormap('hot'); set(gca,'XTick',[], 'YTick', [])

%% visualizing one frame/time instant refocused at different depths
idx_frame = 60;
image_f_stack = squeeze(im_crop(:,:,idx_frame,1:Nstep-1));
im_fm = image_energy(image_f_stack);  % focus measure for the focal stack
[~, idx_debug] = max(im_fm);

figure; montage(norm1(image_f_stack)); colormap('hot'); title('focal stack');
figure; plot(im_fm, 'LineWidth', 5); title('focal measure across stack')

%% As in DfF, filtering the depth map (across time dimension)
idx_depth_filt = medfilt1(idx_depth,9);
figure; subplot(1,2,1); plot(idx_depth_filt,'LineWidth', 5);

depth_vol_flt = depth_estimate(idx_depth_filt);
% l1 trend filtering of the depth map across time
options.niter = 100;
options.method = 'chambolle';
options.lambda = 1e3;
options.verb = 0;
[depth_vol_flt,~] = l1tf(depth_vol_flt.', 400);
subplot(1,2,2); plot(depth_vol_flt - min(depth_vol_flt),'LineWidth', 5);

%% generate 3D video for the results visulaization in Paraview
N_start = 1;
depth_vol_vis = depth_vol_flt(N_start:end);
depth_vol_vis = depth_vol_vis - min(depth_vol_vis);
depth_scale = length(depth_vol_vis)/max(depth_vol_vis);
depth_vol_idx = round(depth_vol_vis*depth_scale);
nx = size(im_crop,1);
DC_vol = 0;
for K = 1: length(depth_vol_vis)
    vis_volume = zeros(nx,nx,length(depth_vol_vis)+20);
    vis_volume(:,:,depth_vol_idx(K)+10) = im_crop(:,:,N_start + K -1,idx_depth(N_start+K-1));
    vis_volume = single(norm1(flip(vis_volume,3)));
    DC_vol = DC_vol + vis_volume;
    vtkwrite([DIR filename '\' filename '_' num2str(K) '.vtk'], 'structured_points', '3D_vol', vis_volume);
end
vtkwrite([DIR filename '\' filename '_DC.vtk'], 'structured_points', '3D_vol', DC_vol);
im_DC_vis = squeeze(sum(DC_vol,3)); 