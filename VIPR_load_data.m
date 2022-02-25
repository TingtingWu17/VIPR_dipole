
%% load data images
load('DH_noise_image_demo.mat');
imageStack = Ix;

XY_data = zeros(size(imageStack,3),2);
Z_data = (0:100:700)*10^-9;

