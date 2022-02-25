
function IS = pre_comp_opt_v2(IS,vec_model_pol)


    lambda = IS.lambda;%wavelength
    M = IS.M;
    pix_size = IS.Cam_psize;
    n1 = IS.n_glass;
    nh = IS.n_glass;
    n2 = IS.n_med;
    NA = IS.NA;
    N = IS.sampling_size;
    zh = 0;

    %calculate both pupil and image plane sampling, 
    %one will affect the other, so make sure not to introduce aliasing

    dx_true = (pix_size/M);%image plane sampling
    dx = n1*dx_true;%due to Abbe sine condition, scale by imaging medium r.i. (see appendix of my journal club)
    %%???????  



    dv = 1/(N*dx);%pupil sampling, related to image plane by FFT
    % recall comb function, 1/dx is the preriod in fourier space(pupil space)
    %%%????

    %define pupil coordinates
    temp=linspace((-1/(2*dx)),(1/(2*dx)),N);
    [eta,xi] = meshgrid(temp);
    
    
    % [eta,xi] = meshgrid(((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(2*N*dx))+(1/(2*dx))),...
    %     ((-1/(2*dx))+(1/(2*N*dx))):dv:(-(1/(N*2*dx))+(1/(2*dx))));

    xBFP = lambda*eta;  % why scaled by wavelength??
    yBFP = lambda*xi;
    [phi,rho] = cart2pol(xBFP,yBFP);
    rho_max = NA/n1;%pupil region of support determined by NA and imaging medium r.i.

    k1 = n1*(2*pi/lambda);
    kh = nh*(2*pi/lambda);
    k2 = n2*(2*pi/lambda);
    %rho(rho >= rho_max) = 0;
    theta1 = asin(rho);%theta in matched medium
    thetah = asin((n1/nh)*sin(theta1));%theta in thin film
    theta2 = asin((n1/n2)*sin(theta1));%theta in mismatched medium
    theta2 = real(theta2)+1i*abs(imag(theta2));


    %%%%%%%%% Start
    %Fresnel coefficients
    tp_2h = 2*n2*cos(theta2)./(n2*cos(thetah) + nh*cos(theta2));
    ts_2h = 2*n2*cos(theta2)./(nh*cos(thetah) + n2*cos(theta2));
    tp_h1 = 2*nh*cos(thetah)./(nh*cos(theta1) + n1*cos(thetah));
    ts_h1 = 2*nh*cos(thetah)./(n1*cos(theta1) + nh*cos(thetah));

    rp_2h = (n2*cos(theta2) - nh*cos(thetah))./(n2*cos(theta2)+ nh*cos(thetah));
    rs_2h = (nh*cos(theta2) - n2*cos(thetah))./(nh*cos(theta2)+ n2*cos(thetah));
    rp_h1 = (nh*cos(thetah) - n1*cos(theta1))./(nh*cos(thetah)+ n1*cos(theta1));
    rs_h1 = (n1*cos(thetah) - nh*cos(theta1))./(n1*cos(thetah)+ nh*cos(theta1));

    %Axelrod's equations for E-fields in back focal plane

   tp = tp_2h.*tp_h1.*exp(1i*kh*cos(thetah)*zh)./(1 + rp_2h.*rp_h1.*exp(2i*kh*zh*cos(thetah))); 
    ts = ts_2h.*ts_h1.*exp(1i*kh*cos(thetah)*zh)./(1 + rs_2h.*rs_h1.*exp(2i*kh*zh*cos(thetah)));

    % Es = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*(muy.*cos(phi) - mux.*sin(phi));
    % Ep = tp.*((n1/n2).*(mux.*cos(phi) + muy.*sin(phi)).*cos(theta1) - muz*sin(theta1).*(n1/n2)^2.*(cos(theta1)./cos(theta2)));
    % Esx - Es contributed by mux

    % still has problem on how to get these to equation????????

    %based on the equation above, seperating the mux,muy,muz compoment from two
    %polirized electric field
    Esx = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*(-sin(phi));
    Esy = ts.*(cos(theta1)./cos(theta2)).*(n1/n2).*cos(phi);
    Epx = tp.*(n1/n2).*cos(phi).*cos(theta1);
    Epy = tp.*(n1/n2).*sin(phi).*cos(theta1);
    Epz = tp.*(-sin(theta1).*(n1/n2)^2.*(cos(theta1)./cos(theta2)));

    % Exx - Ex contributed by mux
    % the first x represents x channel and y channel on the camera, the second
    % x,y,z represents the compoment of mux, muy, muz from the orientation of
    % dipole
    Exx = (1./sqrt(cos(theta1))).*(cos(phi).*Epx - sin(phi).*Esx); %added defocus aberration + depth aberration
    Exy = (1./sqrt(cos(theta1))).*(cos(phi).*Epy - sin(phi).*Esy); %added defocus aberration + depth aberration
    Exz = (1./sqrt(cos(theta1))).*(cos(phi).*Epz); %added defocus aberration + depth aberration
    Eyx = (1./sqrt(cos(theta1))).*(cos(phi).*Esx + sin(phi).*Epx);
    Eyy = (1./sqrt(cos(theta1))).*(cos(phi).*Esy + sin(phi).*Epy);
    Eyz = (1./sqrt(cos(theta1))).*(sin(phi).*Epz);



    % remove the electric component that is outside the accecptant region of
    % objective lens
    Exx(rho >= rho_max) = 0;
    Exy(rho >= rho_max) = 0;
    Exz(rho >= rho_max) = 0;
    Eyx(rho >= rho_max) = 0;
    Eyy(rho >= rho_max) = 0;
    Eyz(rho >= rho_max) = 0;

    %% Oumeng's coordinate flipping
    % coord flipping
    Exx = rot90(Exx);
    Exy = rot90(Exy);
    Exz = rot90(Exz);
    Eyx = fliplr(Eyx);
    Eyy = fliplr(Eyy);
    Eyz = fliplr(Eyz);
    
    %opt_phase_mat(:,:,1) = exp(1i*k1*z*cos(theta1)).*exp(1i*kh*zh*cos(thetah)).*exp(1i*k2*z2*cos(theta2));
    temp2= k1*cos(theta1); % for focus position
    temp1 = k2*cos(theta2);  % for nominal focus plane
    temp1(rho >= rho_max) = 0;
    temp2(rho >= rho_max) = 0;
    
    opt_phase_mat(:,:,1) = temp1;
    opt_phase_mat(:,:,2) = temp2;
    
    IS.opt_phase_mat = opt_phase_mat;
    g_bfpx(:,:,1) = Exx;
    g_bfpx(:,:,2) = Exy;
    g_bfpx(:,:,3) = Exz;
    g_bfpy(:,:,1) = Eyx;
    g_bfpy(:,:,2) = Eyy;
    g_bfpy(:,:,3) = Eyz;
    
    if vec_model_pol == 'x'
        IS.g_bfp = g_bfpx;
    elseif vec_model_pol == 'y'
        IS.g_bfp = g_bfpy;
    end

    %IS.common_term = common_term;
    IS.pz = k1*cos(theta1);
    IS.pz2 = k2*cos(theta2);
    IS.sin_theta = sin(theta1);
    IS.cos_theta = cos(theta1);
    IS.cos_theta2 = cos(theta2);
    IS.XI = xi;
    IS.ETA = eta;


    
end