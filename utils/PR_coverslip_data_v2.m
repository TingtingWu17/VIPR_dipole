function [maskRec,gBlur,Nph,I_mod] = PR_coverslip_data_v2(IS,IMG_T,q,std_stack,gpu_flag,vec_model_flag,cost_function_flag,plot_flag,Alg_flag,est_gBlur_flag,noisy_flag,vec_model_pol,initial_pmask_name)

%% inputs:
% IS : structure with the optical parameters and hyperparameter(defined in main function)
% IMG_T : N images of the z_stack (cropped in a square odd grid) for PR
% q_cord : (x,y,z,NFP) position [um] per image (matrix size NX4)
% gpuflag :  1 - use GPU , 0 - use CPU (10X slower)
% vec_model_flag :   1 - vectorial model, 0 - scalar model
% cost_function_flag :  1 - L1, 2 - L2, 3 - Poiss MLE, 4 - Sum of gaussians MLE
% plot_flag  :  plot while SGD runs, slows down run time X3
% Alg_flag :  1 - ADAM, 2 - Adamax , 3- Nadam, 4 - Nesterov, 5- SVRG, 6- Vanilla SGD
% est_gBlur_flag : flag to estimate the blur in the optimization
% noisy_flag % 0- for dpcesign PSFs, 1 - for PR

%% output
% maskRec - unwraped retrieved phase mask, real pixel size.
% gBlur - updated (if enables) gBlur parameter
% Nph - estimated signal per slice
% I_mod - simulated z-stack, if plot_flag 

%% pre compute all the optical parameters
IS = pre_comp_opt_v2(IS,vec_model_pol);

I_mod=[];
%% needed derivative parameters
step_size = IS.step_size;
k = (2*pi/IS.lambda);
%maskInit = IS.ETA*0;
%load('pixOL_mask.mat');
load(initial_pmask_name);
maskInit =pixOL_mask;
% match the sampling size
maskTemp = zeros(IS.sampling_size,IS.sampling_size);
if IS.sampling_size<size(maskInit,1)
    range = (0:IS.sampling_size-1)-round((IS.sampling_size-1)/2);
    range = range+round(size(maskInit,1)/2);
    maskTemp = maskInit(range,range);
else
    range = (0:size(maskInit,1)-1)-round((size(maskInit,1)-1)/2)+round(IS.sampling_size/2);
    maskTemp(range,range)=maskInit;
end
maskInit = maskTemp;
%maskInit =mask_opt;
%int_cos = IS.int_cos;
z_emit = (IS.z_emit);
opt_phase_mat = (IS.opt_phase_mat);


if vec_model_flag==1
    g_bfp = (IS.g_bfp);
else
    g_bfp = 1;
end
%% move data to GPU if enabled
if gpu_flag == 1
        
    q = gpuArray(q);
    IMG_T = gpuArray(IMG_T);
    std_stack = gpuArray(std_stack);
    step_size = gpuArray(step_size);
    
    k = gpuArray(k);
    
    std_stack = gpuArray(std_stack);
    
    int_cos = gpuArray(int_cos);
    opt_phase_mat = gpuArray(opt_phase_mat);
    
    % medium circ mask
    circ = gpuArray(circ);
    circ_sc = gpuArray(circ_sc);
    
 
end

%% create cost function
% create GPU generator without blur
gb_tmp = IS.gBlur;
IS.gBlur = 0;
IS.FOV_size_crop = IS.FOV_size;

IS.gBlur = gb_tmp;

PSF_mask = @(phase_mask,q) PSF_generator_v2(phase_mask,q,IS);
% define GPU cost function
Cost_fun = @(phase_mask,q,data,std_stack_gpu,Nph_gpu,gBlur_gpu,Nph_opt_flag,cost_altr) CostPR_gpu_v2(phase_mask,IS,q,data,(std_stack_gpu),Nph_gpu,cost_function_flag...
    ,gBlur_gpu ,k, opt_phase_mat ,g_bfp,vec_model_flag,Nph_opt_flag,cost_altr,vec_model_pol,maskInit);

% define GPU init - cost function
Cost_fun_init = @(q,data,std_stack_gpu,Nph_gpu,gBlur_gpu) CostPR_gpu_init_v2(IS,q,data,(std_stack_gpu),Nph_gpu...
     ,k, opt_phase_mat,g_bfp,vec_model_flag,vec_model_pol);
%
% run optimization
[maskRec,gBlur,Nph] = GDM_phasemask_SGDdirect_gpu_v2(IMG_T,Cost_fun,(maskInit),PSF_mask,IS, Alg_flag,plot_flag,est_gBlur_flag,q,std_stack,Cost_fun_init,step_size,gpu_flag,noisy_flag);
Nph = gather(Nph);

%% new generator definition
IS.gBlur = gBlur;
IS.FOV_size = size(IMG_T,1);
PSF_mask = @(phase_mask,q) PSF_generator_v2(phase_mask,q,IS);

if plot_flag==1
    %% plot corr  
    figure(42)
    thr_corr = IS.thr_corr;
    I_mod = [];
    I_dat = [];
    for z_ind = 1:size(q,1)
        tmp = real(PSF_mask(maskRec,q(z_ind,:)));
        I_mod(:,:,z_ind) = gather(tmp);
        %             tmp = tmp(end/2-5+2:end/2+5,end/2-5+2:end/2+5);
        tmp_data = IMG_T(:,:,(z_ind));
        I_dat(:,:,z_ind) = gather(tmp_data);
        % thr mask
        mask = tmp_data>max(tmp_data(:)).*thr_corr;
        tmp = tmp.*mask;
        tmp_data = tmp_data.*mask;
        % normalized corr
        tmp = tmp./norm(tmp(:));
        tmp_data = tmp_data./norm(tmp_data(:));
        corr_coarse(z_ind) = corr2((tmp),(tmp_data));
    end
    plot(q(:,2),corr_coarse)
    title('correlation between rec and data')
    xlabel(' NFP position [um]');
    ylabel('corr2');
    set(gca,'FontSize',10)
    
    %% plot PSF per z position
    imn = @(im) (im - min(im(:)))./(max(im(:)) - min(im(:)));
    for z = 1:size(q,1)
        
        tmp = imn(I_mod(:,:,z));
        tmp = tmp(round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2),round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2));
        tmp_dat = I_dat(:,:,z);
        tmp_dat = tmp_dat(round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2),round(end/2-IS.plotsize/2)+1:round(end/2+IS.plotsize/2));
        
        figure(100)
        subplot(1,2,1);
        imagesc(imn(tmp));
        daspect([1 1 1]);
        title(['model NFP=',num2str(q(z,2))]);
        subplot(1,2,2);
        imagesc(imn(tmp_dat));
        daspect([1 1 1]);
        title(['data NFP=',num2str(q(z,2))]);
        drawnow
        pause(0.1)
        set(gca,'FontSize',10)
    end
    close(100)
end
end
