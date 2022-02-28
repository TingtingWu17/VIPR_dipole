%%
% written by Boris Ferdman
% this functions performs phase retrieval via a direct numerical gradient
% associated paper: VIPR: Vectorial Implementation of Phase Retrieval for
% fast and accurate microscopic pixel-wise pupil estimation 2020

% License to use this code is granted to all interested, as long as the original author is
% referenced as such. The original author maintains the right to be solely associated with this work.
% Copyright by Boris Ferdman: borisferd@gmail.com
%

%clear all;close all;clc
%add current folders file into path
%% data set
VIPR_user_input_dipole;

%%
% don't let plotsize be larger than data
IS.plotsize = min(IS.plotsize,IS.FOV_size);
% define positions per image - [x,y,z,NFP]
q_cord = [z_pos';z_stack_pos']';
IS.psi = q_cord;

%% Gaussian noise estimation per Z
if noisy_flag
    dx=IS.corner_size;
    for j = 1:size(IMG_T,3)
        tmp = [IMG_T(1:dx,1:dx,j);IMG_T(end-dx+1:end-1+1,end-dx+1:end-1+1,j);IMG_T(end-dx+1:end-1+1,1:dx,j);IMG_T(1:dx,end-dx+1:end-1+1,j)];
        %
        mean_stack(j) = mean(tmp(tmp>0));
        % Estimate the std
        std_stack(j) = std(tmp(tmp>0));
        
    end
    IMG_bu = IMG_T;
else
    IMG_bu = IMG_T;
    std_stack = zeros(1,size(q_cord,1));
    mean_stack = zeros(1,size(q_cord,1));
end
%% thr and reduce offset from the images
IMG_T = IMG_bu;
%
for j = 1:size(IMG_T,3)
    tmp = IMG_T(:,:,j)-mean_stack(j);
    % threshold on %
    if IS.I_thr_flag==1
        maskT = tmp>IS.I_thr.*max(tmp(:));
    else
        maskT = tmp>IS.I_thr.*std_stack(j);
    end
    % erode&dilate mask
    se = strel('disk',1,6);
    erodedI = imerode(double(maskT),se);
    se = strel('disk',1,6);
    mask(:,:,j) = imdilate(erodedI,se);
    
    IMG_T(:,:,j) = tmp.*mask(:,:,j);
    
    ff=figure(11)
    subplot(1,3,1)
    imagesc(tmp); colorbar();
    axis  image
    title('input image');
    set(gca,'FontSize',10)
    subplot(1,3,2)
    imagesc(IMG_T(:,:,j)); colorbar();
    axis  image
    title('thresholded image');
    set(gca,'FontSize',10)
    subplot(1,3,3)
    imagesc(IMG_T(:,:,j)-tmp); colorbar();
    axis  image
    title('diff');
    set(gca,'FontSize',10)
    
end
%close(11)


%% flip the image based on the simdipole code
for j=1:size(IMG_T,3)
    tmp = IMG_T(:,:,j);
    if vec_model_pol=='x'
        tmp = flipud(fliplr(tmp)).';
        %tmp = rot90(rot90(tmp,2),3);

    elseif vec_model_pol == 'y'
        %tmp  = fliplr(flipud(tmp));
        tmp  = flipud(tmp);
      
    end
    IMG_T(:,:,j)=tmp;
end
%% Fresnel matching of the grid

mask_size = IS.sampling_size;

%% pad data to match simulation size
tmp_stack = [];
for j = 1:size(IMG_T,3)
    tmp = padarray(IMG_T(:,:,j), floor([mask_size - (size(IMG_T,1)), mask_size - (size(IMG_T,2))]/2),'both');
    if mod(mask_size,2)==0
        tmp = padarray(tmp, floor([1,1]),'pre');
    end
    if mod(size(tmp,1),2)==0
        tmp = padarray(tmp, floor([1,1]),'pre');
    end
    % complete 0 to mean noise
    tmp_stack(:,:,j) = tmp;
    tmp = [];
end
% IMG_T = tmp_stack;

%% retrieve the mask
[maskRec,gB,Nph,I_mod] = PR_coverslip_data_v2(IS,tmp_stack,q_cord,std_stack,gpu_flag,vec_model_flag,cost_function_flag,plot_flag,Alg_flag,est_gBlur_flag,noisy_flag,vec_model_pol,initial_pmask_name);
