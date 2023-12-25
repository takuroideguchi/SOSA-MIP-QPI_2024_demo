clear all
close all

%% get file names
path_sample ='inputs\raw_averaged_images\'; %path to the folder containing the files of the measured images
files_sample = dir([path_sample,'*.csv']); %filenames of the measured images
num_rot = 8; %number of illumination angles measured
num_spect = length(files_sample)/2/num_rot; %number of spectral points measured
path_bg = 'inputs\background_images\'; %path to the foloder containing the measured background images for each illumination angle
files_bg = dir([path_bg,'*.tiff']); %filenames of the measured background images

%% set physical parameters
wav = 532*10^-9; %probe wavelength
NA = 1.2; %NA of the objective lens
pixelsize = 3.45*10^-6/(200/3)/5; %sampling resolution in the sample plane
dim = 1536; %image size
center_x = 650;center_y = 765; %position of the carrier frequency of the interferogram in the Fourier domain
d_pupil_modulation = 27; %diameter of the pupil area affected by the non-sample-specific MIR modulation.
p_pupil_modulation = [1,0]; %position of the modulated pupil area with respect to the objective lens' aperture

%% calculate necessary physical parameters
freq_per_pixel = 1/(pixelsize*dim); %frequency resolution in the Fourier domain (in 1/m)
d_aperture = 2*round((2*NA/wav)/freq_per_pixel/2)+1; %diameter of the objective lens' circular aperture in the Fourier domain (in pixels)

%% create a mask used to crop the effective area of the interferometric term of each interferogram
x = (1:d_aperture) - round(d_aperture/2);
[x,y] = meshgrid(x,x);
mask_aperture = (sqrt(x.^2 + y.^2) <=(d_aperture-1)/2);
mask_pupil_modulation = circshift((sqrt(x.^2 + y.^2) >(d_pupil_modulation-1)/2), p_pupil_modulation);
mask = mask_aperture.*mask_pupil_modulation;

%% read background images
for k = 1:num_rot
    temp_img  = double(imread(strcat(path_bg,files_bg(k).name))); %background hologram without any sample in the field of view
    temp_holo_bg_stack(:,:,k) = temp_img / mean2(temp_img);
end

%% perform synthetic-aperture reconstruction
for k = 1:num_spect
    for l = 1:num_rot
        %read holograms
        temp_holo_1 = transpose(csvread(strcat(path_sample,files_sample(2*num_rot*(k-1)+(l-1)*2+1).name))); %hologram of the sample corresponding to either MIR ON or OFF state
        temp_holo_2 = transpose(csvread(strcat(path_sample,files_sample(2*num_rot*(k-1)+(l-1)*2+2).name))); %hologram of the sample corresponding to the other state (i.e., MIR OFF or ON)
        temp_holo_bg = temp_holo_bg_stack(:,:,l);

        % normalize the brightness of the holgoram
        temp_holo_1 = temp_holo_1/mean2(temp_holo_1);
        temp_holo_2 = temp_holo_2/mean2(temp_holo_2);
    
        %calculate frequency spectra of the holograms
        temp_freq_1 = fftshift(fft2(temp_holo_1));
        temp_freq_2 = fftshift(fft2(temp_holo_2));
        temp_freq_bg = fftshift(fft2(temp_holo_bg));

        % crop the interferometric term
        temp_freq_1 = temp_freq_1(center_y-(d_aperture-1)/2:center_y+(d_aperture-1)/2,center_x-(d_aperture-1)/2:center_x+(d_aperture-1)/2).*mask_aperture;
        temp_freq_2 = temp_freq_2(center_y-(d_aperture-1)/2:center_y+(d_aperture-1)/2,center_x-(d_aperture-1)/2:center_x+(d_aperture-1)/2).*mask_aperture;
        temp_freq_bg = temp_freq_bg(center_y-(d_aperture-1)/2:center_y+(d_aperture-1)/2,center_x-(d_aperture-1)/2:center_x+(d_aperture-1)/2).*mask_aperture;

        % zeropad the frequency spectra so that the reconstructed images are upsampled by a factor of 2.
        temp_freq_1 = padarray(temp_freq_1,[round(d_aperture/2) round(d_aperture/2)],0);
        temp_freq_2 = padarray(temp_freq_2,[round(d_aperture/2) round(d_aperture/2)],0);
        temp_freq_bg = padarray(temp_freq_bg,[round(d_aperture/2) round(d_aperture/2)],0);
        temp_reconst_1_raw = ifft2(ifftshift(temp_freq_1))./ifft2(ifftshift(temp_freq_bg));
        temp_reconst_2_raw = ifft2(ifftshift(temp_freq_2))./ifft2(ifftshift(temp_freq_bg));       

        % calcluate the MIR ON-OFF differential image.
        % This image contains artifact associated with the non-sample-specific MIR modulation of the pupil function
        temp_reconst_MIRdiff_raw = temp_reconst_1_raw./temp_reconst_2_raw;

        % calculate the frequency spectra of the three reconstructions
        temp_freq_reconst_1_raw = fftshift(fft2(temp_reconst_1_raw));
        temp_freq_reconst_2_raw = fftshift(fft2(temp_reconst_2_raw));
        temp_freq_reconst_MIRdiff_raw = fftshift(fft2(temp_reconst_MIRdiff_raw));

        % estimate the shift vector of the oblique illumination in the Fourier
        % domain. Based on this, calculate the Fourier mask that is used to
        % crop the effective area of the Fourier spectrum of hte recontructed
        % complex field that is unaffected by the non-sample-specific MIR
        % modulation. Apply these masks.
        [temp,idx] = max(abs(temp_freq_1));
        [temp,idx_x] = max(temp);
        idx_y = idx(idx_x);
        temp_mask = circshift(padarray(mask,[round(d_aperture/2) round(d_aperture/2)],0),[d_aperture-idx_y,d_aperture-idx_x]);
        temp_freq_reconst_1 = temp_freq_reconst_1_raw.*temp_mask;
        temp_freq_reconst_2 = temp_freq_reconst_2_raw.*temp_mask;
        temp_freq_reconst_MIRdiff = temp_freq_reconst_MIRdiff_raw.*temp_mask;
        
        % calculate the MIR ON-OFF differential image after eliminating the
        % non-sample-specific MIR modulation of the pupil function
        temp_reconst_MIRdiff =  ifft2(ifftshift(temp_freq_reconst_MIRdiff));

        % check if the central part of the differential image shows increase or decrease in
        % phase. Depending on the result, determine if temp_freq_reconst_1 or
        % temp_freq_reconst_2 is MIR ON or OFF state.
        temp_reconst_MIRdiff = temp_reconst_MIRdiff/mean2(temp_reconst_MIRdiff(5:2*d_aperture+1-5,5:2*d_aperture+1-5));
        if mean2(angle(temp_reconst_MIRdiff(60:100,60:80)))<0
            temp_reconst_MIRON = ifft2(ifftshift(temp_freq_reconst_1));
            temp_reconst_MIROFF = ifft2(ifftshift(temp_freq_reconst_2));
        else
            temp_reconst_MIRON = ifft2(ifftshift(temp_freq_reconst_2));
            temp_reconst_MIROFF = ifft2(ifftshift(temp_freq_reconst_1));
        end 
       
        % remove the reconstruction artifact that may appear at the edges
        % of the images
        x = -d_aperture:d_aperture;
        [x,y] = meshgrid(x,x);
        edge_size = 2;
        temp_edge_to_average = (max(abs(x),abs(y)) == d_aperture+1-6);
        temp_edge_to_remove = (max(abs(x),abs(y)) >= d_aperture+1-edge_size);
        temp_reconst_MIROFF(temp_edge_to_remove)=sum(sum(temp_reconst_MIROFF(temp_edge_to_average)))/sum(sum(temp_edge_to_average));
        temp_reconst_MIRON(temp_edge_to_remove)=sum(sum(temp_reconst_MIRON(temp_edge_to_average)))/sum(sum(temp_edge_to_average));

        % calculate the three reconstructions
        temp_reconst_MIROFFON = temp_reconst_MIROFF./temp_reconst_MIRON;
        temp_reconst_MIROFF = temp_reconst_MIROFF/temp_reconst_MIROFF(1,1);
        temp_reconst_MIRON = temp_reconst_MIRON/temp_reconst_MIRON(1,1);

        % calculate the frequency spectra of the reconstructions
        temp_freq_reconst_MIROFFON = fftshift(fft2(temp_reconst_MIROFFON)).*temp_mask;
        temp_freq_reconst_MIROFF = fftshift(fft2(temp_reconst_MIROFF)).*temp_mask;
        temp_freq_reconst_MIRON = fftshift(fft2(temp_reconst_MIRON)).*temp_mask;

        % collect temporary outputs into variables
        freq_MIROFFON_stack(:,:,k,l) = temp_freq_reconst_MIROFFON;
        freq_MIROFF_stack(:,:,k,l) = temp_freq_reconst_MIROFF;
        freq_MIRON_stack(:,:,k,l) = temp_freq_reconst_MIRON;
        mask_stack(:,:,k,l) = temp_mask;  

        figure(1);
        subplot(1,3,1);imagesc(log(abs(temp_freq_1)));daspect([1 1 1])
        subplot(1,3,2);imagesc(log(abs(temp_freq_2)));daspect([1 1 1])
        subplot(1,3,3);imagesc(log(abs(temp_freq_bg)));daspect([1 1 1])
        figure(2);
        subplot(1,3,1);imagesc(angle(temp_reconst_MIRdiff_raw));daspect([1 1 1]);
        subplot(1,3,2);imagesc(angle(temp_reconst_MIRdiff));daspect([1 1 1]);
        subplot(1,3,3);imagesc(angle(temp_reconst_MIROFFON));daspect([1 1 1]);
        figure(3);
        subplot(1,4,1);imagesc(temp_mask);daspect([1 1 1])
        subplot(1,4,2);imagesc(log(abs(temp_freq_reconst_MIROFFON)));daspect([1 1 1])
        subplot(1,4,3);imagesc(log(abs(temp_freq_reconst_MIROFF)));daspect([1 1 1])
        subplot(1,4,4);imagesc(log(abs(temp_freq_reconst_MIRON)));daspect([1 1 1])
        figure(4);imagesc(angle(temp_reconst_MIROFF));daspect([1 1 1])
        drawnow
        l
    end

    % calculate the synthetic aperture from all the angular measurements.
    % Also upsample the synthetic spectrum by a factor of 4.
    temp_freq = sum(squeeze(freq_MIROFFON_stack(:,:,k,:)),3)./(sum(squeeze(mask_stack(:,:,k,:)),3)+((sum(squeeze(mask_stack(:,:,k,:)),3)==0)));
    temp_freq_upsampled = padarray(temp_freq,4*[d_aperture, d_aperture],0);
    freq_synth_MIROFFON(:,:,k) = temp_freq;
    freq_synth_MIROFFON_upsampled(:,:,k) = temp_freq_upsampled;
    reconst_synth_MIROFFON(:,:,k) = ifft2(ifftshift(temp_freq));
    reconst_synth_MIROFFON_upsampled(:,:,k) = ifft2(ifftshift(temp_freq_upsampled));

    temp_freq = sum(squeeze(freq_MIROFF_stack(:,:,k,:)),3)./(sum(squeeze(mask_stack(:,:,k,:)),3)+((sum(squeeze(mask_stack(:,:,k,:)),3)==0)));
    temp_freq_upsampled = padarray(temp_freq,4*[d_aperture, d_aperture],0);
    freq_synth_MIROFF(:,:,k) = temp_freq;
    freq_synth_MIROFF_upsampled(:,:,k) = temp_freq_upsampled;
    reconst_synth_MIROFF(:,:,k) = ifft2(ifftshift(temp_freq));
    reconst_synth_MIROFF_upsampled(:,:,k) = ifft2(ifftshift(temp_freq_upsampled));

    temp_freq = sum(squeeze(freq_MIRON_stack(:,:,k,:)),3)./(sum(squeeze(mask_stack(:,:,k,:)),3)+((sum(squeeze(mask_stack(:,:,k,:)),3)==0)));
    temp_freq_upsampled = padarray(temp_freq,4*[d_aperture, d_aperture],0);
    freq_synth_MIRON(:,:,k) = temp_freq;
    freq_synth_MIRON_upsampled(:,:,k) = temp_freq_upsampled;
    reconst_synth_MIRON(:,:,k) = ifft2(ifftshift(temp_freq));
    reconst_synth_MIRON_upsampled(:,:,k) = ifft2(ifftshift(temp_freq_upsampled));

    figure(5);
    subplot(2,3,1);imagesc(log(abs((freq_synth_MIROFFON(:,:,k)))));daspect([1 1 1])
    subplot(2,3,4);imagesc(angle(reconst_synth_MIROFFON_upsampled(:,:,k)));daspect([1 1 1])
    subplot(2,3,2);imagesc(log(abs((freq_synth_MIROFF(:,:,k)))));daspect([1 1 1])
    subplot(2,3,5);imagesc(angle(reconst_synth_MIROFF_upsampled(:,:,k)));daspect([1 1 1])
    subplot(2,3,3);imagesc(log(abs((freq_synth_MIRON(:,:,k)))));daspect([1 1 1])
    subplot(2,3,6);imagesc(angle(reconst_synth_MIRON_upsampled(:,:,k)));daspect([1 1 1])
    drawnow;
    k
end
figure(6);imagesc(angle(reconst_synth_MIROFF_upsampled(:,:,1)));daspect([1 1 1]);colormap('gray');caxis([-0.3,0.8]);colorbar
saveas(gcf,'outputs/OPD_MIROFF_2920.fig')
figure(7);imagesc(angle(reconst_synth_MIROFFON_upsampled(:,:,1)));daspect([1 1 1]);colorbar
saveas(gcf,'outputs/MIP_2920.fig')
figure(8);imagesc(angle(reconst_synth_MIROFF_upsampled(:,:,2)));daspect([1 1 1]);colormap('gray');caxis([-0.3,0.8]);colorbar
saveas(gcf,'outputs/OPD_MIROFF_3150.fig')
figure(9);imagesc(angle(reconst_synth_MIROFFON_upsampled(:,:,2)));daspect([1 1 1]);colorbar
saveas(gcf,'outputs/MIP_3150.fig')