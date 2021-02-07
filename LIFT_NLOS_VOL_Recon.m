% ************* This script perform LIFT reconstruciton for NLOS *******************
% ****** It performs the reconstruction in a volumetric (all frames) manner ********
% ****** the (x,y,t) datacube is not unwarped here ***********
close all;
clear all;
clc;
top_DIR = '.\NLOS\' 
filename = 'LetterVI'; % LetterN
lambda_TV = 1e-1;    % the hyperparameter for the FISTA reconstruction
IMG_INV = false;  % use image relay where image is inverted
dt = 4.6e-12;    % temporal resolution: 4.6 ps or 10 ps for TwoCircles and Circle_Square

% ******
% 1: Load the raw data with system calibration results
% ******
mat_content = load([top_DIR '\raw_data_' filename '.mat']);
Calib_Res = mat_content.Calib_Res;
streak_tform = Calib_Res.streak_tform;
image_t = mat_content.image_t;
% ******
% 2: Load image raw data and correct for the streak distortion 
% ******
figure(1); subplot(2,2,1); imagesc(image_t); colormap('hot'); title('raw streak image');
image_t = imwarp(image_t, streak_tform, 'OutputView', imref2d(size(image_t)));
subplot(2,2,2); imagesc(norm1(image_t)); title('rectified streak image');
drawnow();

% *******
% 3a: setup solver and refocusing parameters
% *******
u = (1:7)-4; s = 0; % shearing factor for refocusing
% sub_img_cnt = round(Calib_Res.sub_img_cnt + u.*s); % normal refoucsing procedure
% disabling the dpeth-fof-field version: supplementary note 3.3
sub_img_cnt = round(Calib_Res.sub_img_cnt + u.*(s.*tan(deg2rad(90-Calib_Res.Angle))).');

options = []; options.INVERT = IMG_INV; options.CROP = true;
options.Deconv = false; options.USE_TV = false;
options.Refocus = true;
options.sub_img_cnt = round(sub_img_cnt); 
options.Normalize = true;

% 3b: call the solver
tic
im_crop = fx_LIFT_ReconVOL(Calib_Res, image_t, lambda_TV, options);
toc
im_crop = gather(im_crop);

subplot(2,2,3); imagesc(sum(im_crop,3)); title('DC image');

%%
% % *********
% % 4: unwarp the NLOS data
% % *********
nlos_data = nlos_unwarp(im_crop, Calib_Res, dt);
DeltaT = 3e8*dt;
laser_pos = Calib_Res.laser_pos;  % laser pos on the wall, relative to global coo.
grid_pos = Calib_Res.grid_pos; 
ReceiveDelay = Calib_Res.ReceiveDelay;
% save the results for NLOS reconstrcution
save([top_DIR '\NLOS_' filename '.mat'], 'nlos_data', ...
    'DeltaT','grid_pos','laser_pos','ReceiveDelay');