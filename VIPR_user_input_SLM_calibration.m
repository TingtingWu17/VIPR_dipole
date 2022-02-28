%% user defined flags
gpu_flag = 0; % 1 - use GPU, 0 - on CPU
vec_model_flag = 1; % 1 - vectorial model, 0 - scalar model
cost_function_flag = 4; % optimization cost 1 - L1, 2 - L2, 3 - Poiss MLE, 4 - Sum of gaussians MLE 5 - Sum of gaussians MLE+phase mask difference
plot_flag = 1; % plot while SGD runs, slows down ~ X4
Alg_flag = 1  ; % gradient method : 1 - ADAM, 2 - Adamax , 3- Nadam, 4 - Nesterov ,5- Vanilla SGD
vec_model_pol = 'x' ; %'x' or 'y' for having a  polarizer, 'b' for full vectorial
noisy_flag = 1; % 0- for design PSFs, 1 - for PR;
est_gBlur_flag = 1; % 1- estimate gBlur after 1/3 of the iterations
crop_flag = 1; % 1 = point choice, 2= CoG of max projection , 3- no cropping
%initial_pmask_name = 'pixOL_mask_com.mat';
pmasm_name_bmp = 'checkerboard_000.bmp';
pmask_image = double(imread(pmasm_name_bmp));
pmask_image = pmask_img/256*2*pi;

%% define optical parameters
IS.Signal=1; % for generator - keep 1
IS.bg=0; % for generator - keep 0
%IS.tophoton = 0.29;
IS.tophoton = 0.66;

%optical parameters
IS.M=111.1111; % objective magnification
IS.NA=1.40; % objective NA
IS.lambda=676e-9; % wavelength [um]
IS.Cam_psize = 6.5*10^-6; % pixel size of camera [m]
IS.gBlur=0.75; %initial guess of the blurring
IS.n_glass=1.518; % RI of immersion oil
IS.n_med=1.518; % RI of sample medium
%IS.f_4f = 15e4; % focal length of 4-f  system (if 2D imaging- use tube f)

IS.SAF_flag = 1; % include SAF or not (1=include - recommended)
IS.Iris_pNA = 1; % optional iris to limit BFP, in range [0,1] where 1 is the NA
% emitter size
IS.z_emit = 50*10^-9; % emitter radius [m]
IS.bfp_radius = 80;
IS.sampling_size = round(IS.lambda*IS.bfp_radius/(IS.Cam_psize/IS.M)/IS.NA);
if rem(IS.sampling_size,2)==0
    IS.sampling_size=IS.sampling_size+1;
end
% polarization of dipole (0,0,0) is for freely rotating
IS.p_vec = [0,0,0]; % incoherent


%% load data
% two example datas: the first one is captured for all NFP focused at
% beads; the second one is capture for NFP axially scanning the beads
if exampleData == 1
beadsDatafile  = [pwd '\TW_pixOL_com_beads_data\20220214_combine_beads_data_152_to_data_164.mat']; % note the beads data: y channel is same as images captured on the camera, that is there is no need to flip the y channel image
load(beadsDatafile);
IS.FOV_size = size(beads_img,1); % size of ROI used
if vec_model_pol=='x'
    SM1 = beads_img(:,[1:IS.FOV_size],:)*IS.tophoton;
elseif  vec_model_pol=='y'
    SM1 = beads_img(:,[1:IS.FOV_size]+IS.FOV_size,:)*IS.tophoton;
end
IMG_T  = SM1;
xy = zeros(size(IMG_T,3),2);
z_stack_pos = zeros(size(SM1,3),1)-90*10^-9; %position of NFP
z_pos = zeros(size(SM1,3),1)+IS.z_emit;  % position of the emitter

elseif exmapleData ==2
beadsDatafile  = [pwd '\TW_pixOL_com_beads_data\20210602_data12_beads_y_channel_unflipped.mat'];

load(beadsDatafile);
IS.FOV_size = size(beads_img,1); % size of ROI used
if vec_model_pol=='x'
    SM1 = beads_img(:,[1:IS.FOV_size],:)*IS.tophoton;
elseif  vec_model_pol=='y'
    SM1 = beads_img(:,[1:IS.FOV_size]+IS.FOV_size,:)*IS.tophoton;
end
IMG_T  = SM1;
xy = zeros(size(IMG_T,3),2);
z_stack_pos = [(-90-[-700:50:-50,50:50:700,0])*10^-9];
z_stack_pos = repmat(z_stack_pos,11,1);
z_stack_pos = reshape(z_stack_pos,[],1);
z_pos = zeros(size(SM1,3),1)+IS.z_emit;  % position of the emitter
end
%% optimization parameters
% the parameters in this section define the optimization proccess

% pre-proc parameters
IS.I_thr_flag = 2; % 1- thr above IS.thr*max(I) per image, else - thr above IS.thr*background_std
IS.I_thr = 1; %  threshold parameter
IS.corner_size = 5; % size of corners to  estimate noise [pixels]

% hyper-params
IS.SGDiter = 250; %  how  many iterations to SGD
IS.step_size = 3e-1; % step  size (try 3e-1 for ADAM and 3e-8 for SGD)
IS.point_num = 9; % size of mini-batch per SGD iteration

% additional options
IS.gBlur_cost = 4; % cost to estimate gBlur if est_gBlur_flag=1, (1-4 same as cost_function_flag , 5 - by corr2)
% option not to use SGD
IS.last_iter = 0; % how many iterations  to run not with SGD (at end of optimization)
IS.last_iter_flag = 3; % 1 - contuine SGD, 2 - global gradient, 3- batch the lowest correlation points, 
%                        4- adaptive sampling with side info on corr
IS.thr_corr = 0.01; % threshold for correlation calc (used if last_iter_flag = 3)
IS.upsample_fact = 1; % how much to upsample the data (usually leave at 1)
IS.update_Signal = 1; % 1 - update signal at second half of iterations 
IS.count_Nph_epochs=3; % how  many times to update signal(if IS.update_signal ==1) - too many might over-fit
IS.Photobleach_cost_mult = 0; % add to cost the SNR consideration
% plot sizes for PSF
IS.plotsize = 99 ; % size of psf plots [pixels]
