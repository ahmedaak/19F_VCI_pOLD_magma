% This script analyzes FAIR_EPI (inversion recovery T1 mapping) for the 19F
% (pOLD) in VCI study 

% Requirements: FSL,

% Written by Ahmed Khalil 15.9.2017 (+ using tools by Dr. Philipp Boehm-Sturm)

% add paths
addpath(genpath('/home/khalila/VCI/Ahmed/Ahmed_FittingScripts'));

clear all
datadir = '/home/khalila/VCI/Ahmed/19F_VCI_Study/19F_VCI_POLD/GC_optimize/';
cd(datadir)
mouse_folder=dir('M*');

%mouse_folder = uigetdir('\19F_VCI_pOLD\BASELINE\dat\GC_optimize\', 'Select folder with FAIR-EPI data')
TIs=[25 400 800 1200 1600 2000]; % define inversion times used in FAIR-EPI

for dirnum=1:length(mouse_folder);
cd(mouse_folder(dirnum).name)

fair_epi_list = dir('*FAIR_EPI_T1map*');
if exist('output')
    system('rm -r output');
end 

if length(fair_epi_list) == 0
    disp('There is no FAIR-EPI data in this folder')
elseif length(fair_epi_list) > 1
    disp('There are multiple FAIR-EPI datasets in this folder')
else 
    disp('Analyzing a single FAIR-EPI dataset from this folder')
end

for f = 1:2; % loop through fitting modes (1=all voxels, 2=only those with SNR>5)
for h=1:2; % loop through SNR mask types (from T2 vs FAIR-EPI)
    if h==1;
        snr_mask_name = dir('SNR_T2_*')
        snr_mask_name = snr_mask_name.name;
        snr_type = 'T2mask';
    elseif h==2;
        snr_mask_name = dir('SNR_FAIR_*')
        snr_mask_name = snr_mask_name.name;
        snr_type = 'FAIRmask';
    end
% loop through FAIR-EPI datasets
params_mean_all = [];
errs_mean_all = [];
for i=1:length(fair_epi_list)

      fair_epi_name = fair_epi_list(i).name;  
        k = findstr(fair_epi_name, '30');
        l = findstr(fair_epi_name, '60');
        m = findstr(fair_epi_name, '100');
        
        if ~isempty(k);
        o2_level = '30';
    elseif ~isempty(l);
        o2_level = '60';
    elseif ~isempty(m);
        o2_level = '100';
        else
            disp('Please make sure the FAIR-EPI file is named properly, i.e. with the O2 level (30, 60, or 100)')
    end
        system(['fslsplit ' fair_epi_name ' -z']); % split data slicewise
        fair_epi_4D_name = strrep(fair_epi_name, '.nii', '_4D.nii'); % new name for merged (4D) data
        system(['fslmerge -t ' fair_epi_4D_name ' vol0*']);  % concatenate data timewise
        system('rm vol0*') % delete split files
        
        fair_epi_4D_sm = strrep(fair_epi_4D_name, '.nii', '_sm.nii'); % define name of smoothed data
        system(['fslmaths ' fair_epi_4D_name ' -kernel gauss 0.6 -fmean ' fair_epi_4D_sm]) % smooth data
        snr_mask_name_bin = strrep(snr_mask_name, '.nii', '_bin.nii');
        system(['fslmaths ' snr_mask_name ' -bin ' snr_mask_name_bin]); % binarize SNR mask
        
        system(['fslcpgeom ' fair_epi_4D_sm ' ' snr_mask_name_bin]); % fix orientation issue (R/L flipping) - WILL MAKE MASK 4D - SEE NEXT STEP
        system(['fslroi ' snr_mask_name_bin ' ' snr_mask_name_bin ' 0 1']); % take only first "volume" of mask
        
        system(['fslmaths ' snr_mask_name_bin ' -mul FAIR_crude_brainmask.nii.gz ' snr_mask_name_bin]) % "brain extract" SNR mask
        mean_signal_file = ['mean_signal_' o2_level '.txt'];
        system(['fslmeants -i ' fair_epi_4D_sm ' -m ' snr_mask_name_bin ' -o ' mean_signal_file]); % extract raw signal values from SNR mask (average)
      

        
        if f == 2; % perform SNR thresholding   
              
        % remove all voxels outside SNR mask
        fair_epi_4D_thresh = strrep(fair_epi_4D_sm, '.nii', [snr_type, '_thr.nii']); % define thresholded FAIR-EPI file name

        system(['fslmaths ' fair_epi_4D_sm ' -mul ' snr_mask_name_bin ' ' fair_epi_4D_thresh]); % perform thresholding
        system(['fslcpgeom ' fair_epi_4D_sm ' ' fair_epi_4D_thresh]); % fix orientation issue (R/L flipping)
        outputdir = [snr_type '_thr'];
        
        elseif f == 1; % Do not perform SNR thresholding
        fair_epi_4D_thresh = fair_epi_4D_sm;
        outputdir = [snr_type];

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
    save_nii(params_nii,[o2_level, 'params.nii']);
    err_nii=make_nii(err);
    save_nii(err_nii,[o2_level, 'errs.nii']);
    
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
    system(['fslsplit ' [o2_level, 'params.nii'] ' -z']); % split data slicewise
    system(['mv vol0000.nii.gz ' o2_level '_T1.nii.gz; mv vol0001.nii.gz ' o2_level '_S0.nii.gz; mv vol0002.nii.gz ' o2_level '_beta.nii.gz']);
    % fix orientation issue with SNR mask BEFORE extracting values (fixes
    % R/L flipping issue)
    system(['fslcpgeom 30_T1.nii.gz ' snr_mask_name_bin]);
    % extract and save parameter values
    system(['fslmeants -i ' o2_level '_T1.nii.gz -m ' snr_mask_name_bin ' --showall --transpose -o ' o2_level '_T1.txt' ])
    system(['fslmeants -i ' o2_level '_S0.nii.gz -m ' snr_mask_name_bin ' --showall --transpose -o ' o2_level '_S0.txt' ])
    system(['fslmeants -i ' o2_level '_beta.nii.gz -m ' snr_mask_name_bin ' --showall --transpose -o ' o2_level '_beta.txt' ])

end

if ~exist(['output/'])
system('mkdir output');
end
system(['mkdir output/' outputdir]);
system(['mv *_sm.* *_T1.* *_beta.* *_S0.* *4D.nii.gz *thr.nii *thr.nii.gz *bin.nii.gz *.txt *params.nii *errs.nii *3D.nii output/' outputdir]);
end
end
cd(datadir)
end
    