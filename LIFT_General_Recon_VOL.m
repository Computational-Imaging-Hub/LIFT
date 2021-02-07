% ************* This script perform LIFT reconstruciton for general Time domain measurments *******************
% ****** It performs the reconstruction in a volumetric (all frames) manner ********
close all;
clear all;
clc;

top_DIR = '.\LIF\DOF_Ver\'; 
filename = 'LinearFiber_vertical';
lambda_TV = 7e-1;   % the TV regularization parameter
IMG_INV = false;

% ******
% 1: Load the calibration data resultsfigure
% ******
Calib_Res = load([top_DIR,'Calib_Res.mat']);
Calib_Res = Calib_Res.Calib_Res;
streak_tform = Calib_Res.streak_tform;

% reconFile_HIS = [top_DIR filename '\' filename '.his'];  %
% image_t = readHIStoOne(reconFile_HIS,'SUM');
% image_t = image_t - mean( mean(image_t(900:1000,:)) );
% save([top_DIR filename '\' filename '.mat']);

% ******
% 2: Load image data and extract the slice useful for image reconstruction
% ******
image_t = load([top_DIR filename '\' filename '.mat']);
image_t = image_t.image_t;

figure; subplot(2,2,1); imagesc(image_t);
axis square;  title('Streak image'); colormap('hot');

% ******
%  3a: Correct for the streak distortion 
% ******
image_t = imwarp(image_t, streak_tform, 'OutputView', imref2d(size(image_t)));

subplot(2,2,2); imagesc(norm1(image_t)); title('rectified streak image');
drawnow();

% crop the signal to reduce reconstruction time: time slices full of 0s are not useful
image_t = image_t(150:1:450,:);

%% 3: setup solver
options = []; options.INVERT = IMG_INV;  options.CROP = true; options.Deconv = true;
options.Normalize = false; options.USE_TV = false;
options.Refocus = true; % refocus flag
options.sub_img_cnt = round(Calib_Res.cntx_depth(:,5)).';
tic
im_crop = fx_LIFT_ReconVOL(Calib_Res, image_t, lambda_TV, options);
toc
%%
im_DC = sum(im_crop,3);  % the time-integrated image
subplot(2,2,3); imagesc(im_DC); axis square; title('reconstructed DC image');
