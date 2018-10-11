% This script analyzes FAIR_EPI (inversion recovery T1 mapping) for the 19F
% (pOLD) in VCI study - in vitro calibration data

% Requirements: FSL,

% Written by Ahmed Khalil 14.05.2018 (+ using tools by Dr. Philipp Boehm-Sturm)

% add paths
addpath(genpath('/home/khalila/VCI/Ahmed/19F_VCI_Study/19F_VCI_POLD'));

clear all
datadir = '/home/khalila/VCI/Ahmed/19F_VCI_Study/19F_VCI_invitro_calibration_April2018/dat_analyzed_May18/19F_VCI_POLD/in_vitro_calibration/dat/';
cd(datadir)
mouse_folder=dir('20180425_vciPOLDcalib_*');

TIs=[25 400 800 1200 1600 2000]; % define inversion times used in FAIR-EPI

for dirnum=1:length(mouse_folder); 
    cd(mouse_folder(dirnum).name)

    fair_epi_list = dir('*FAIR_EPI_T1map*');
    if exist('output')
        system('rm -r output');
    end 

        for f = 1:2; % loop through fitting modes (1=all voxels, 2=only those with SNR>5)
        % loop through FAIR-EPI datasets
        params_mean_all = [];
        errs_mean_all = [];
        for i=1:length(fair_epi_list)

              fair_epi_name = fair_epi_list(i).name;  
              snr_mask_name = 'SNR_FAIREPI.nii';

        % PERFORM PREPROCESSING
                system(['fslsplit ' fair_epi_name ' -z']); % split data slicewise
                fair_epi_4D_name = strrep(fair_epi_name, '.nii', '_4D.nii'); % new name for merged (4D) data
                system(['fslmerge -t ' fair_epi_4D_name ' vol0*']);  % concatenate data timewise
                system('rm vol0*') % delete split files

                fair_epi_4D_sm = strrep(fair_epi_4D_name, '.nii', '_sm.nii'); % define name of smoothed data
                system(['fslmaths ' fair_epi_4D_name ' -kernel gauss 0.3 -fmean ' fair_epi_4D_sm]) % smooth data
                snr_mask_name_bin = strrep(snr_mask_name, '.nii', '_bin.nii');
                system(['fslmaths ' snr_mask_name ' -bin ' snr_mask_name_bin]); % binarize SNR mask

                system(['fslcpgeom ' fair_epi_4D_sm ' ' snr_mask_name_bin]); % fix orientation issue (R/L flipping) - WILL MAKE MASK 4D - SEE NEXT STEP
                system(['fslcpgeom ' fair_epi_4D_sm '  FAIR_crude_brainmask.nii.gz']); % fix orientation issue (R/L flipping)
                system(['fslmaths FAIR_crude_brainmask.nii.gz -bin  FAIR_crude_brainmask.nii.gz']); % binarize brain mask 

                system(['fslroi ' snr_mask_name_bin ' ' snr_mask_name_bin ' 0 1']); % take only first "volume" of mask
                system(['fslroi  FAIR_crude_brainmask.nii.gz FAIR_crude_brainmask.nii.gz 0 1']); % take only first "volume" of mask

                system(['fslmaths ' snr_mask_name_bin ' -mul FAIR_crude_brainmask.nii.gz ' snr_mask_name_bin]) % "brain extract" SNR mask
                mean_signal_file = ['mean_signal.txt'];
                system(['fslmeants -i ' fair_epi_4D_sm ' -m ' snr_mask_name_bin ' -o ' mean_signal_file]); % extract raw signal values from SNR mask (average)


                if f == 2; % Perform SNR thresholding   

                % remove all voxels outside SNR mask
                fair_epi_4D_thresh = strrep(fair_epi_4D_sm, '.nii', '_thr.nii'); % define thresholded FAIR-EPI file name

                system(['fslmaths ' fair_epi_4D_sm ' -mul ' snr_mask_name_bin ' ' fair_epi_4D_thresh]); % perform thresholding
                system(['fslcpgeom ' fair_epi_4D_sm ' ' fair_epi_4D_thresh]); % fix orientation issue (R/L flipping)
                outputdir = 'SNR_FAIREPI';

                elseif f == 1; % Do NOT perform SNR thresholding
                fair_epi_4D_thresh = fair_epi_4D_sm;
                outputdir = 'noSNR';

                else disp('Not a valid fitting mode');
                end


                % split data timewise
                system(['fslsplit ' fair_epi_4D_thresh ' -t']);
                % concatenate data slicewise
                fair_epi_3D_thresh = strrep(fair_epi_4D_thresh, '4D', '3D');
                system(['fslmerge -z ' fair_epi_3D_thresh ' vol0*']);
                system('rm vol0*') % delete split files
                system(['fslchfiletype NIFTI ' fair_epi_3D_thresh]) % convert to noncompressed NIFTI


                %% Perform fitting VOXELWISE  
            % load in NIFTI data
            rawdata_nii=load_nii(fair_epi_3D_thresh);
            rawdata_array=rawdata_nii.img;    

            [params,err]=fit_nlinfit(rawdata_array,TIs,'T1_invRec',0);
            %plotTC(rawdata_array,1,'T1_invRec',TIs);
            params_nii=make_nii(params);
            save_nii(params_nii,['params.nii']);
            err_nii=make_nii(err);
            save_nii(err_nii,['errs.nii']);

                %% Perform fitting MEAN SIGNAL IN SNR MASK
                % load mean signal
            meansig = dlmread(mean_signal_file);
                % reshape into 1x1x6 matrix
            meansig_rs = reshape(meansig, [1 1 6]);
                % fit
            [params_mean, err_mean] = fit_nlinfit(meansig_rs, TIs, 'T1_invRec',0);
                % save parameters
                params_mean_all = [params_mean_all; params_mean(1)];
                errs_mean_all = [errs_mean_all; err_mean(1)];
             csvwrite('t1_mean.txt',params_mean_all);

            % split parameter files into single files (slicewise) and rename
            system(['fslsplit ' ['params.nii'] ' -z']); % split data slicewise
            system(['mv vol0000.nii.gz ' 'T1.nii.gz; mv vol0001.nii.gz ' 'S0.nii.gz; mv vol0002.nii.gz ' 'beta.nii.gz']);
            % fix orientation issue with SNR mask BEFORE extracting values (fixes
            % R/L flipping issue)
            system(['fslcpgeom 30_T1.nii.gz ' snr_mask_name_bin]);
            % extract and save parameter values
            system(['fslmeants -i ' 'T1.nii.gz -m ' snr_mask_name_bin ' --showall --transpose -o ' 'T1.txt' ])
            system(['fslmeants -i ' 'S0.nii.gz -m ' snr_mask_name_bin ' --showall --transpose -o ' 'S0.txt' ])
            system(['fslmeants -i ' 'beta.nii.gz -m ' snr_mask_name_bin ' --showall --transpose -o ' 'beta.txt' ])

        end

        if ~exist(['output/'])
        system('mkdir output');
        end
        system(['mkdir output/' outputdir]);
        system(['mv *_sm.* *T1.* *beta.* *S0.* *4D.nii.gz *thr.nii *thr.nii.gz *bin.nii.gz *.txt *params.nii *errs.nii *3D.nii output/' outputdir]);
        end
cd(datadir)
end
    