function PSF = PSF_generator_v2(pmask,q,IS)
    z0 = q(1);
    NFP = q(2);
    opt_phase_mat = IS.opt_phase_mat;
    BFP_Phase = (z0.*opt_phase_mat(:,:,1)+NFP*opt_phase_mat(:,:,2));
    BFP_Phase = exp(1i.*(BFP_Phase+pmask));
    g_bfp = IS.g_bfp;
    
    pmask = exp(1i*pmask);
    imgEx = fftshift(fft2(g_bfp(:,:,1).*BFP_Phase));
    imgEy = fftshift(fft2(g_bfp(:,:,2).*BFP_Phase));
    imgEz = fftshift(fft2(g_bfp(:,:,3).*BFP_Phase));


%     % coord flipping
%     imgExx = fliplr(flipud(imgExx'));
%     imgExy = fliplr(flipud(imgExy'));
%     imgExz = fliplr(flipud(imgExz'));
%     imgEyx = flipud(imgEyx);
%     imgEyy = flipud(imgEyy);
%     imgEyz = flipud(imgEyz);

    % euqation from backer's paper Eq.22
    basisImage(:,:,1) = abs(imgEx).^2;
    basisImage(:,:,2) = abs(imgEy).^2;
    basisImage(:,:,3) = abs(imgEz).^2;
    basisImage(:,:,4) = 2*real(conj(imgEx).*imgEy);
    basisImage(:,:,5) = 2*real(conj(imgEx).*imgEz);
    basisImage(:,:,6) = 2*real(conj(imgEy).*imgEz);  
    %the results are same with the above equation





    PSF = 1/3*basisImage(:,:,1)+1/3*basisImage(:,:,2)+1/3*basisImage(:,:,3);
        

    PSF = PSF./sum(sum(PSF));

end