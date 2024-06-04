function rcbFLOW()
%
% rcbFLOW: rCBF fitting, lineup, and observation wizard
%   
%     Passning av perfusionsdata enligt Kober at al. NMR Biomed
%     21(8):781-92. Calculates regional cerebral blood flow (rCBF) from 
%     flow-sensitive alternating inversion recovery echo planar imaging (FAIR-EPI) 
%     acquisitions over a user-defined region of interest (ROI) mask in one brain slice 
%     per batch data set and automatically registers the resultant masked rCBF image
%     to the average over a user-defined stack of corresponding anatomical
%     reference image slices. Updated for MATLAB 2023b. Parallel Processing
%     Toolbox required. 
%
%     Add to MATLAB path before calling rcbFLOW: 
% 
%          rcbFLOW.m
%          rcbFLOW_invrecexp.m
%          rcbFLOW_DMTI2.m
%          rcbFLOW_DMTIb.m
%
%     Inputs (by message prompt): 
%           
%           1   batch_processing (string; '0' or '1'): '0' processes a single
%               dataset; '1' enables batch analysis of a list of inputs defined as
%               in the example CSV below, each additionally defined in the case of
%               batch_processing = 0
%           2   auto_masking for FAIR-EPI (string; '0', '1', or '2'): '0' requires the user to manually 
%               select the analysis ROI on a sample FAIR-EPI image, while
%               '1' uses a multi-step segmentation process to do so without
%               manual user input. '2' enables loading of previously traced
%               masks saved as dcms. 
%           3   auto_masking for RARE (string; '0', '1', or '2'): '0' requires the user to manually 
%               select the analysis ROI on a sample RARE image, while
%               '1' uses a multi-step segmentation process to do so without
%               manual user input. '2' enables loading of previously traced
%               masks saved as dcms. 
%           4   anti_aliasing (string; '0' or '1'): '0' will not implement
%               anti-aliasing correction on rCBF images, while '1' will. 
%
%     Outputs (per every dataset analyzed, by file export to single new directory): 
%           
%           1   *_rCBF.dcm (DCM): Masked rCBF slice calculated from FAIR-EPI
%           2   *_T1c.dcm (DCM): Masked T1c slice calculated from FAIR-EPI
%           3   *_RARE_uncorrected.dcm (DCM): Averaged anatomical reference slice
%           4   *_RARE_flipped_resized.dcm (DCM): Flipped and resized
%               anatomical reference slice
%           5   *_RARE_flipped_resized_registered.dcm (DCM): Flipped, resized, and registered 
%               anatomical reference slice
%           6   FAIR_EPI_mask_origsize_*.dcm (DCM): Mask for FAIR-EPI (rCBF
%               and T1c) slice
%           7   RARE_mask_origsize_*.dcm (DCM): Mask for RARE reference slice
%           8   RARE_mask_resized_*.dcm (DCM): Mask for RARE reference
%               slice resized for overlay on rCBF and T1c slices 
%           9  *_Fig1 (PNG): Figure 1. Difference between selective and global
%               inversion conditions for each inversion time (TI) 
%           10  *_Fig2 (PNG): Figure 2. FAIR-EPI ROI mask (either manual, automatic, or loaded) for rCBF
%               calculation 
%           11  *_Fig3 (PNG): Figure 3. RARE ROI mask (either manual, automatic, or loaded)
%               for registration
%           12  *_Fig4 (PNG): Figure 4. Masked calculated R1c image
%           13  *_Fig5 (PNG): Figure 5. Masked calculated rCBF image 
%           14  *_Fig6 (PNG): Figure 6. Masked calculated T1c image
%           15  *_Fig7 (PNG): Figure 7. Averaged anatomical reference image 
%           16  *_Fig8 (PNG): Figure 8. Flipped anatomical reference image 
%           17  *_Fig9 (PNG): Figure 9. Masked calculated rCBF image overlaid on averaged anatomical reference image 
%           18  *_Fig10 (PNG): Figure 10. Masked rCBF image overlaid on flipped, registered, and resized averaged anatomical reference image 
%           
%     Example batch input CSV: 
% 
%           FAIR_directory_fullpath,RARE_directory_fullpath,first_RARE_slice_filename, (optional: FAIR_mask_file_fullpath, RARE_mask_file_fullpath)
%           C:\Users\kswanberg\CBFDir1,C:\Users\kswanberg\RAREDir1,MRIm12.dcm,C:\Users\kswanberg\FAIRmasks\1.dcm,C:\Users\kswanberg\RAREmasks\1.dcm,
%           C:\Users\kswanberg\CBFDir2,C:\Users\kswanberg\RAREDir2,MRIm14.dcm,C:\Users\kswanberg\FAIRmasks\2.dcm,C:\Users\kswanberg\RAREmasks\2.dcm,
%           C:\Users\kswanberg\CBFDir3,C:\Users\kswanberg\RAREDir3,MRIm10.dcm,C:\Users\kswanberg\FAIRmasks\3.dcm,C:\Users\kswanberg\RAREmasks\3.dcm,
%           ...
%           C:\Users\kswanberg\CBFDir75,C:\Users\kswanberg\RAREDir75,MRIm12.dcm,C:\Users\kswanberg\FAIRmasks\75.dcm,C:\Users\kswanberg\RAREmasks\75.dcm,
% 
%     Contributors: 
% 
%     Michael Gottschalk, 202303, conception, image processing, voxelwise curve fitting
%     Kelley Swanberg, 202312-2406, flex I/O, auto image seg/reg, batching, docs
% 
%     Copyright (C) 2023 Michael Gottschalk and Kelley M. Swanberg
% 

%% STEP 1: Define input parameters for file I/O and rCBF fitting 
close all
clear all

% Define currently hard-coded input parameters  
t = [0.03,0.5,1,3,5,9.7]'; % Here we are using FAIR-EPI with 6 TI 
num_ti = length(t); 
tsize = 2*num_ti;
num_slicesR = 3; % Number of RARE slices to average into a reference slice
R1blood = 0.64;
lambda = 90*60; % ml/100g/min Kober, multiply by 60 otherwise per second 
alfa = 0.96;
tau = 2;
H = [];

% Create directory for output files 
slicedirroot = pwd; 
datetime_var = strrep(strrep(strrep(string(datetime), ':', '_'), ' ', '_'), '-', '_');
slicedirname = sprintf("%s\\%s_perfusionRefOutput", slicedirroot, datetime_var);

% Prompt user to process single case or batch 
prompt = {'Will you be processing a batch of cases or a single case today? [Batch=1; Single=0]'};
dlgtitle = 'Input';
dims = [1 32];
definput = {'0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
batch_processing = cell2mat(answer(1));

% Check for a valid batch input 
if ~(strcmp(batch_processing,'0') || strcmp(batch_processing,'1'))
    batch_processing = '0';
    message = sprintf('The choices were [Batch=1; Single=0]. \nSingle=%s will now be used by default.', batch_processing);
    uiwait(warndlg(message, 'Input Warning', 'modal'));
end

% Optional automatic rCBF map masking 
prompt = {'How would you like to mask FAIR-EPI slices for RCBF map calculation? [Load=2; Auto=1; Manual=0]'};
dlgtitle = 'Input';
dims = [1 32];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
auto_masking = cell2mat(answer(1));

% Check for a valid rCBF map masking input 
if ~(strcmp(auto_masking,'0') || strcmp(auto_masking,'1') || strcmp(auto_masking,'2'))
    auto_masking = '1';
    message = sprintf('The choices were [Load Masks=2; Automatic Masking=1; Manual Masking=0]. \nAutomatic Masking=%s will now be used by default.', auto_masking);
    uiwait(warndlg(message, 'Input Warning', 'modal'));
end

% Optional automatic RARE masking 
prompt = {'How would you like to mask RARE slices for anatomical referencing? [Load=2; Auto=1; Manual=0]'};
dlgtitle = 'Input';
dims = [1 32];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
auto_maskingR = cell2mat(answer(1));

% Check for a valid rCBF map masking input 
if ~(strcmp(auto_maskingR,'0') || strcmp(auto_maskingR,'1') || strcmp(auto_maskingR,'2'))
    auto_maskingR = '1';
    message = sprintf('The choices were [Load Masks=2; Automatic Masking=1; Manual Masking=0]. \nAutomatic Masking=%s will now be used by default.', auto_maskingR);
    uiwait(warndlg(message, 'Input Warning', 'modal'));
end

% Ask user whether to implement optional anti-aliasing correction 
prompt = {'Would you like to correct RCBF map for aliasing (recommended)? [Yes=1; No=0]:'};
dlgtitle = 'Input';
dims = [1 32];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
anti_aliasing = cell2mat(answer(1));

% Check for a valid anti-aliasing correction input 
if ~(strcmp(anti_aliasing,'0') || strcmp(anti_aliasing,'1'))
    anti_aliasing = '1';
    message = sprintf('The choices were [Anti-aliasing=1; No anti-aliasing=0]. \nAnti-aliasing=%s will now be used by default.', anti_aliasing);
    uiwait(warndlg(message, 'Input Warning', 'modal'));
end

% Get file inputs for batch or single case 
if batch_processing == '1' % Batch processing 

    % Prompt user for batch input CSV
    uiwait(msgbox('Please select your batch-processing input CSV (Columns: FAIR_directory_fullpath, RARE_directory_fullpath, first_RARE_slice_filename, FAIR_EPI_mask_fullpath (optional), RARE_mask_fullpath (optional))!'))
    [input_CSV_path,~] = uigetfile({'*.csv'}, 'Select your input CSV.', slicedirroot,'MultiSelect','off');
    input_CSV = readtable(input_CSV_path, 'Delimiter',',','NumHeaderLines',0);
    num_datasets = height(input_CSV); 

    % Check that all folders and files from input CSV actually exist and will not cause problems later  
    for dataset_ii = 1:num_datasets

        for inputCSVcol_kk = 1:3

            if inputCSVcol_kk<3
                inputPath = string(input_CSV{dataset_ii, inputCSVcol_kk});

                if ~exist(inputPath, 'dir')
                    error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath)
                end

            else 
                inputPathPath = string(input_CSV{dataset_ii, inputCSVcol_kk-1});
                inputPathFile = string(input_CSV{dataset_ii, inputCSVcol_kk}); 
                inputPath = strcat(inputPathPath, '\', inputPathFile); 

                if ~exist(inputPath, 'file')
                    error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath); 
                end

            end
            
        end

    end 

     % Check mask inputs only if load mask option selected 
    if auto_masking == '2'
        if auto_maskingR == '2' % FAIR-EPI mask path in column 4 and RARE mask path in column 5 
            for dataset_ii = 1:num_datasets
                inputPath = string(input_CSV{dataset_ii, 4});

                if ~exist(inputPath, 'file')
                    error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath)
                end

                inputPath = string(input_CSV{dataset_ii, 5});

                if ~exist(inputPath, 'file')
                    error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath)
                end

            end 
        else % Only FAIR-EPI mask in column 4
            for dataset_ii = 1:num_datasets
                  inputPath = string(input_CSV{dataset_ii, 4});

                if ~exist(inputPath, 'file')
                    error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath)
                end

            end 
        end
    else 
        if auto_maskingR == '2' % Only RARE mask path in column 4 
            for dataset_ii = 1:num_datasets
                inputPath = string(input_CSV{dataset_ii, 4});

                if ~exist(inputPath, 'file')
                    error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath)
                end

            end 
        end
    end


    disp('Initial check of input CSV successful!')

else % Single-case processing 

    num_datasets = 1; 
    % Prompt user for FAIR-EPI acquisitions to process 
    FAIREPI_input_prompt = sprintf('Please select all %d DICOMs from the full FAIR-EPI inversion-recovery series in the same slice!', tsize); 
    uiwait(msgbox(FAIREPI_input_prompt))
    [filename,~] = imgetfile('InitialPath',slicedirroot,'MultiSelect',true);

    % Ensure number of selected files equals twice number of expected
    % inversion times 
    if length(filename) ~= tsize 
        error('You have selected %d FAIR-EPI DICOMs, but we are expecting %d FAIR-EPI DICOMs instead.', length(filename), tsize); 
    end

    % Check that all selected files are dicoms 
    for file_ii = 1:length(filename)
        if ~isdicom(cell2mat(filename(file_ii)))
            error('File %d that you have selected is not a DICOM.', file_ii)
        end
    end

    % Import RARE data according to user input  
    RARE_input_prompt = sprintf('Please select the %d RARE DICOM template files corresponding to the FAIR-EPI slice!', num_slicesR); 
    uiwait(msgbox(RARE_input_prompt))
    [filenameR,~] = imgetfile('InitialPath',slicedirroot,'MultiSelect',true);

    % Ensure number of selected files equals predefined input for number
    % expected
    if length(filenameR) ~= num_slicesR 
        error('You have selected %d RARE DICOMs, but we are expecting %d RARE DICOMs instead.', length(filenameR), num_slicesR)
    end
      
    % Check that all selected files are DICOMs 
    for file_ii = 1:length(filenameR)
        if ~isdicom(cell2mat(filenameR(file_ii)))
            error('File %d that you have selected is not a DICOM.', file_ii)
        end
    end

    % Import masks according to user input 
    if auto_masking == '2' 

        Mask_input_prompt = sprintf('Please select the DICOM mask for the FAIR-EPI slice!'); 
        uiwait(msgbox(Mask_input_prompt))
        [filename_Mask,~] = imgetfile('InitialPath',slicedirroot,'MultiSelect',true);
        filename_Mask = filename_Mask{1};

        % Check that selected files is a DICOM
        if ~isdicom(filename_Mask)
            error('The file that you have selected is not a DICOM.')
        end

    end

    if auto_maskingR == '2'

        MaskR_input_prompt = sprintf('Please select the DICOM mask for the RARE slice!'); 
        uiwait(msgbox(MaskR_input_prompt))
        [filename_MaskR,~] = imgetfile('InitialPath',slicedirroot,'MultiSelect',true);
        filename_MaskR = filename_MaskR{1}

        % Check that selected files is a DICOM
        if ~isdicom(filename_MaskR)
            error('The file that you have selected is not a DICOM.')
        end

    end

end

%% xxxxxxxxxxxxxxxxxxx THE BIG PROCESSING LOOP xxxxxxxxxxxxxxxxxxx %%

for dataset_ii = 1:num_datasets

    %% STEP 2: Load and preprocess FAIR-EPI acquisitions for rCBF fitting
    if batch_processing == '1'
        % Define FAIR-EPI files in this data set 
        dir_fairepi = string(input_CSV{dataset_ii, 1}); 
        filename_struct = dir(fullfile(dir_fairepi, '*.dcm'));
        filename_path = string({filename_struct(:).folder}); 
        filename_files = string({filename_struct(:).name}); 
        filename = (strcat(filename_path', '\', filename_files'))';

        % Define three RARE files for this data set 
        dir_RARE = string(input_CSV{dataset_ii, 2}); 
        filenameR_all = dir(fullfile(dir_RARE, '*.dcm'));
        RARE_firstfile = string(input_CSV{dataset_ii, 3}); 
        RARE_firstfile_index = find(contains(string({filenameR_all(:).name}'),RARE_firstfile));
        RARE_list = string({filenameR_all(RARE_firstfile_index:RARE_firstfile_index+2).name}'); 
        filenameR_path = string({filenameR_all(1:3).folder}); 
        filenameR = (strcat(filenameR_path', '\', RARE_list))'; 

        % Define mask files in this data set 
        if auto_masking == '2' 
            filename_Mask = string(input_CSV{dataset_ii, 4}); % FAIR-EPI mask path in Column 4
            if auto_maskingR == '2' % RARE mask path in Column 5
                filename_MaskR = string(input_CSV{dataset_ii, 5}); 
            end
        else 
            if auto_maskingR == '2' % RARE mask path in Column 4 
                filename_MaskR = string(input_CSV{dataset_ii, 4}); 
            end
        end 
    end 

    % Prepare filenames for saving 
    [a,~,~]= fileparts((string(filename(1))));
    pathparts = strsplit(a,filesep);
    filename_export = sprintf('%s_%d', pathparts(end), dataset_ii); 
    
    % If export directory does not already exist, make it so 
    if not(isfolder(slicedirname))
        mkdir(slicedirname) 
    end 

    % Cycle through the FAIR-EPI image files provided by user and acquire
    % rescale values by slope and intercept (required for Bruker) 
    for ii = 1:num_ti
        % Load appropriate file as well as pair into memory 
        pair_index = 2*ii-1;
        pathperf1 = cell2mat(filename(pair_index));
        pathperf2 = cell2mat(filename(pair_index+1));
    
        % Load header and image information 
        q1 = dicominfo(pathperf1);
        q2 = dicominfo(pathperf2); 
    
        image_original1 = dicomread(pathperf1);
        image_original2 = dicomread(pathperf2);
    
        % When B does not exist initialize matrix according to DICOM header
        % information 
        if ~exist('B', 'var')
            rn = q1.Height; %100 % read direction FAIR-EPI 
            pn = q1.Width; %120 % phase-encode direction FAIR-EPI
            B = zeros(rn,pn,tsize);
            sequencename = q1.SequenceName; 
            if strcmp(sequencename, 'Bruker:FAIR_EPI')
                sequencekey = 'normalFAIR'; 
            else 
                sequencekey = 'presatFAIR'; 
            end
        end

        % Put individual slice information into array of structs and matrix 
        H = [H q1 q2]; 
        B(:,:,pair_index) = (q1.RescaleSlope)*image_original1 + q1.RescaleIntercept;
        B(:,:,pair_index+1) = (q2.RescaleSlope)*image_original2 + q2.RescaleIntercept;
        
        % This is where presatFAIR (t1a FAIR) and normal FAIR (Bruker FAIR)
        % diverge 
        if strcmp(sequencekey, 'normalFAIR')
            if pair_index < 7
                B(:,:,pair_index) = -B(:,:,pair_index); 
                B(:,:,pair_index+1) = -B(:,:,pair_index + 1); 
            end
        else
            if pair_index < 5
                B(:,:,pair_index) = -B(:,:,pair_index); 
                B(:,:,pair_index+1) = -B(:,:,pair_index + 1); 
            end
        end
    
        % Alternating slices are expected to be inverts of each other %%?
        %if rem(i,2)~=0
            %B(:,:,i) = -B(:,:,i); 
        %end
    end
    
    % Initialize matrix for delta B 
    DB = zeros(rn,pn,num_ti);
    
    % Subtract non-selective inversions from selective inversions for each
    % TI because selective should have more signal due to exchange from
    % greater number of unsaturated / uninverted spins 
    figure(1);
    for ii = 1:num_ti
        pair_index = 2*ii-1; 
        DB(:,:,ii)= double(B(:,:,pair_index))-double(B(:,:,pair_index + 1));  
        subplot(3,2,ii);
        imagesc(DB(:,:,ii));
    end

    sgtitle({'Figure 1. FAIR-EPI selective minus global inversion', 'conditions by inversion time (TI)'})
    name_fig1 = sprintf("%s\\%s_Fig1.png", slicedirname, filename_export); 
    saveas(gcf,name_fig1)
    
    %% STEP 3: Mask FAIR-EPI brain ROI to facilitate computational efficiency for fitting as well as reasonable registration to RARE reference 
    % Create a freehand mask of the brain slice 
    Borig = 0.5*(B(:,:,11)+B(:,:,12));
    
    % % Draw the freehand slice mask 
    if auto_masking == '0'
    
        % Display brain slice for freehand masking 
        figure(2); imagesc(B(:,:,12))
        colormap gray
        daspect([1 1 1])
        ax = gca;
        c = ax.Visible;
        %ax.Visible = 'off';
    
        % Draw the freehand slice mask 
        h = drawfreehand;
        BW = createMask(h);
    
        title({'Figure 2. Manually segmented FAIR-EPI ROI', 'for rCBF calculation'})
        name_fig2 = sprintf("%s\\%s_Fig2.png", slicedirname, filename_export); 
        saveas(gcf,name_fig2)
    
        % Walk through freehand mask one pixel at a time and convert image to a
        % binary mask matrix 
        Mask = zeros(size(B(:,:,12)));
        for ii =1:rn
            for jj =1:pn
                if  BW(ii,jj) < 1
                    Mask(ii,jj) = 0;
                else
                    Mask(ii,jj) = 1;
                end          
            end
        end
    
    elseif auto_masking == '2'  % Load mask from DICOM 
    
        Mask = double(dicomread(filename_Mask));

        % Display brain slice with loaded mask 
        figure(2); imagesc(Mask.*B(:,:,12))
        colormap gray
        daspect([1 1 1])
        ax = gca;
        c = ax.Visible;
        %ax.Visible = 'off';
    
        title({'Figure 2. Loaded mask applied to FAIR-EPI ROI', 'for rCBF calculation'})
        name_fig2 = sprintf("%s\\%s_Fig2.png", slicedirname, filename_export); 
        saveas(gcf,name_fig2)

    else % Automatic masking 

        level_step1 = multithresh(B(:,:,12), 6); 
        h_step1 = imquantize(B(:,:,12), level_step1); 
        BW = imfill(h_step1, "holes"); 
    
        % Walk through freehand mask one pixel at a time and convert image to a
        % binary mask matrix 
        Mask = zeros(size(B(:,:,12)));
        for ii =1:rn
            for jj =1:pn
                if  BW(ii,jj) > 2
                    Mask(ii,jj) = 1;
                else
                    Mask(ii,jj) = 0;
                end          
            end
        end
    
        % Display brain slice with automated mask 
        figure(2); imagesc(Mask.*B(:,:,12))
        colormap gray
        daspect([1 1 1])
        ax = gca;
        c = ax.Visible;
        %ax.Visible = 'off';
    
        title({'Figure 2. Auto-segmented FAIR-EPI ROI', 'for rCBF calculation'})
        name_fig2 = sprintf("%s\\%s_Fig2.png", slicedirname, filename_export); 
        saveas(gcf,name_fig2)

    end
    
    %% STEP 4: Load and preprocess three RARE slices to serve as an averaged reference for rCBF slice
    % Figure out what new image size should be based on RARE reference
    % images and initialize BR matrix according to DICOM header information 
    pathperfR = cell2mat(filenameR(1));
    qR = dicominfo(pathperfR); % Acquire dicominfo 
    rnR = qR.Width; %read direction RARE
    pnR = qR.Height; %phase-encode direction RARE
    
    BR = zeros(pnR,rnR,num_slicesR);
    
    % Cycle through the first three slices added by the user 
    for ii = 1:num_slicesR 
        pathperfR = cell2mat(filenameR(ii));
        qR = dicominfo(pathperfR); % Acquire dicominfo 
        imageR = dicomread(pathperfR); 
        BR(:,:,ii) = (qR.RescaleSlope) .*double(imageR) + qR.RescaleIntercept;
    end
    
    % Sum the three slices together 
    ImageR = BR(:,:,1)+BR(:,:,2)+BR(:,:,3);
    
    if auto_maskingR == '0'
        % Display brain slice for freehand masking 
        figure(3); imagesc(ImageR)
        colormap gray
        daspect([1 1 1])
        ax = gca;
        c = ax.Visible;
        %ax.Visible = 'off';
        
        % Draw the freehand slice mask 
        h = drawfreehand;
        BWR = createMask(h);
        
        title({'Figure 3. Manually segmented RARE ROI', 'for registration'})
        name_fig3 = sprintf("%s\\%s_Fig3.png", slicedirname, filename_export); 
        saveas(gcf,name_fig3)

        % Walk through freehand mask one pixel at a time and convert image to a
        % binary mask matrix 
        MaskR = zeros(size(ImageR));
        for ii =1:pnR
            for jj =1:rnR
                if  BWR(ii,jj) < 1
                    MaskR(ii,jj) = 0;
                else
                    MaskR(ii,jj) = 1;
            end          
            end
        end

    elseif auto_maskingR == '2' % Load mask from DICOM  
        MaskR = double(dicomread(filename_MaskR));

        % Display brain slice with loaded mask 
        figure(3); imagesc(MaskR.*ImageR)
        colormap gray
        daspect([1 1 1])
        ax = gca;
        c = ax.Visible;
    
        title({'Figure 3. Loaded mask applied to RARE ROI', 'for registration'})
        name_fig3 = sprintf("%s\\%s_Fig3.png", slicedirname, filename_export); 
        saveas(gcf,name_fig3)

    else % Automatic masking 
        
        level_step1R = multithresh(ImageR, 6); 
        h_step1R = imquantize(ImageR, level_step1R); 
        BWR = imfill(h_step1R, "holes"); 
    
        % Walk through freehand mask one pixel at a time and convert image to a
        % binary mask matrix 
        MaskR = zeros(size(ImageR));
        for ii =1:pnR
            for jj =1:rnR
                if  BWR(ii,jj) > 2
                    MaskR(ii,jj) = 1;
                else
                    MaskR(ii,jj) = 0;
                end          
            end
        end
    
        % Display brain slice with automated mask 
        figure(3); imagesc(MaskR.*ImageR)
        colormap gray
        daspect([1 1 1])
        ax = gca;
        c = ax.Visible;
    
        title({'Figure 3. Auto-segmented RARE ROI', 'for registration'})
        name_fig3 = sprintf("%s\\%s_Fig3.png", slicedirname, filename_export); 
        saveas(gcf,name_fig3); 

    end
    
    % Apply mask to RARE 
    ImageR = MaskR.*ImageR; 

    %% STEP 5: Fit FAIR-EPI images to solve for M_0, a_0, and especially R1_app
    % Fit data with selective instead of global inversion
    R1c = zeros(rn,pn);
    INVc = R1c;
    Mzeroc = R1c;
    
    if strcmp(sequencekey, 'normalFAIR')

        M_ss = zeros(1,num_ti);

        for ii = 1:rn
            parfor jj = 1:pn   
                if Mask(ii,jj) > 0   

                M_ss_worker = zeros(1,num_ti);

                for kk = 1:num_ti
                    pair_index = 2*kk-1; 
                    M_ss_worker(kk) = B(ii,jj,pair_index);
                end
                
                M_ss = M_ss_worker;
                    y = M_ss'./10000; 
                    ft = fittype('rcbFLOW_invrecexp(x,a,b,c)'); % This is the same 
                    options = fitoptions(ft);
                    options.StartPoint = [max(y), 2*max(y), 1];
                    [fitobject,gof] = fit(t,y,ft,options);
                    R1c(ii,jj) = fitobject.c;
                    INVc(ii,jj) = fitobject.b/(2*fitobject.a);
                    Mzeroc(ii,jj) = fitobject.a;

                end
            end
            disp(ii)
        end

    else

        B_norm = B./B(:,:,12);
        
        % Walk through images in B and correct for infinite values because
        % otherwise our fit will not converge 
        for ii = 1:12
            for jj = 1:pn
                for kk = 1:rn
                    if B_norm(kk,jj,ii) == Inf 
                        B_norm(kk,jj,ii) = 0;
                    elseif isnan(B_norm(kk,jj,ii))
                        B_norm(kk,jj,ii) = 0;
                    else
                    end
                end
            end
        end
        
        % Fit the single inversion-recovery function M_ss = M_0*(1-2*a_0*e^(-TI/R1_app)) to solve for M_0, a_0, and R1_app  
        M_ss = zeros(1,num_ti);
        for ii = 1:rn
        
            parfor jj = 1:pn % All voxels in same phase-encode step are processed in parallel 
        
                if Mask(ii,jj) > 0       
        
                    % Need to initialize another matrix to enable parallel
                    % processing (cannot just use M_ss)
                    M_ss_worker = zeros(1,num_ti);
    
                    % Find M_ss values corresponding to this voxel from each TI 
                    for kk = 1:num_ti
        
                        pair_index = 2*kk-1; 
                        M_ss_worker(kk) = B_norm(ii,jj,pair_index);
        
                    end
        
                    M_ss = M_ss_worker
                    y = M_ss_worker'; 
                    ft = fittype('rcbFLOW_invrecexp(x,a,b,c)');
                    options = fitoptions(ft);
                    options.StartPoint = [1, 2, 0.5];
                    options.MaxIter = 1000; 
                    options.Lower = [-1e16, -1e16, 0]; 
                    [fitobject,gof] = fit(t,y,ft,options);
                    R1c(ii,jj) = fitobject.c;
                    INVc(ii,jj) = fitobject.b/(2*fitobject.a);
                    Mzeroc(ii,jj) = fitobject.a;
        
                end
                
            end
        
            disp(ii)
        
        end
    end

    figure(4);
    imagesc(R1c,[0 0.8])
    
    title({'Figure 4. Local longitudinal relaxation rate', 'R1 estimated from FAIR-EPI'})
    name_fig4 = sprintf("%s\\%s_Fig4.png", slicedirname, filename_export); 
    saveas(gcf,name_fig4)
    
    %% STEP 6: Fit R1c images solved in previous step for final rCBF images
    % Set variable inputs for M
    E = DB;

    % Prepare DeltaM and R1fact variables for linear fit
    if strcmp(sequencekey, 'normalFAIR')
            scalefact = 10000;

            parfor ii = 1:num_ti
            DeltaM(:,:,ii) = lambda.*E(:,:,ii)./(2.*INVc.*scalefact.*Mzeroc);
            R1fact(:,:,ii)= (exp(-t(ii).*R1c)-exp(-t(ii).*R1blood))./(R1blood-R1c);  
            end
            
            RCBF = zeros(rn,pn);
            
            %linear fit including offset to allow for differences in starting point
            %(also present in Kober et al. see fig 4)
            for ii = 1:rn
                parfor jj = 1:pn
                    if Mask(ii,jj) > 0

                    wt = zeros(1,num_ti);
                    tit = zeros(1,num_ti);
                    for kk = 1:num_ti
                        wt(kk) = DeltaM(ii,jj,kk);
                        tit(kk) = R1fact(ii,jj,kk);
                    end
                    w = wt;
                    ti = tit;

                        ft = fittype('rcbFLOW_DMTIb(x,a,b)');
                        options = fitoptions(ft);
                        options.StartPoint = [100 1];
                        options.Lower = [1 -1e16]; 
                        fitobject = fit(ti',w',ft,options);
                        RCBF(ii,jj) = fitobject.a;  
                    end
                end
                disp(ii)
            end
    else

        Mzeroc = Borig.*Mzeroc;
        parfor ii = 1:num_ti
            %DeltaM(:,:,ii) = lambda.*E(:,:,ii)./(2.*alfa.*Mzeroc.*(1-exp(-tau.*R1blood))); % Follows the paper, except why the extra (1-exp(-tau.*R1blood)))
            DeltaM(:,:,ii) = lambda.*E(:,:,ii)./(2.*INVc.*Mzeroc.*(1-exp(-tau.*R1blood)));
            R1fact(:,:,ii)= (exp(-t(ii).*R1c)-exp(-t(ii).*R1blood))./(R1blood-R1c);
        end

        % Use DeltaM and R1fact expressions to fit rCBF pixel by pixel according to
        % equation in Kober et al.
        RCBF = zeros(rn,pn);
        for ii = 1:rn
            
            % Parallel computing ftw 
            parfor jj = 1:pn
                
                if Mask(ii,jj) > 0
        
                    wt = zeros(1,num_ti);
                    tit = zeros(1,num_ti);
        
                    for kk = 1:num_ti
                        wt(kk) = DeltaM(ii,jj,kk);
                        tit(kk) = R1fact(ii,jj,kk);
                    end
        
                    w = wt;
                    ti = tit;
                    % ft = fittype('rcbFLOW_DMTI2(x,a)');
                    % options = fitoptions(ft);
                    % options.StartPoint = 100;
                    ft = fittype('rcbFLOW_DMTIb(x,a,b)');
                    options = fitoptions(ft);
                    options.StartPoint = [100 1];
                    options.Lower = [1 -1e16]; 
                    fitobject = fit(ti',w',ft,options);
                    RCBF(ii,jj) = fitobject.a; 
                end          
            end
            disp(ii)
        end
    end
    
    Mask_output = Mask; 

    %% STEP 7: Anti-aliasing correction for final R1/T1 and rCBF images; output as DCM 
    %shift image down 1 pixel
    for i = 2:rn
        RCBF1(i,:) = RCBF(i-1,:);
        R1c1(i,:) = R1c(i-1,:);
        Mask_output(i,:) = Mask(i-1,:); 
    end
        RCBF1(1,:) = RCBF(rn,:);
        R1c1(1,:) = R1c(rn,:);
        Mask_output(1,:) = Mask(rn,:); 

    if anti_aliasing =='1'
        q1.ImagePositionPatient(1) = q1.ImagePositionPatient(1)+q1.PixelSpacing(1)*double(q1.AcquisitionMatrix(1))*(double(q1.AcquisitionMatrix(1))/double(q1.Width)-1);
    end
    
    figure(5);
    imagesc(RCBF1,[0 0.65*max(RCBF1(:))])
    
    title({'Figure 5. Final rCBF image', 'estimated from FAIR-EPI'})
    name_fig5 = sprintf("%s\\%s_Fig5.png", slicedirname, filename_export); 
    saveas(gcf,name_fig5)

    % Save uncorrected rCBF DICOM 
    filename_dcm = sprintf("%s\\%s_rCBF.dcm", slicedirname, filename_export); 
    dicomwrite(int16(RCBF1),filename_dcm,q1);

    %% STEP 8: Preregister RARE reference to upsampled T1w image from rCBF fit 
    
    % Register RARE to T1 image calculated from rCBF  
    % Invert R1 image
    T1c1=1./R1c1; 
    % Walk through T1c1 and correct for infinite values 
    for jj = 1:rn
        for kk = 1:pn
            if T1c1(jj,kk) == Inf 
                T1c1(jj,kk) = 0;
            elseif isnan(T1c1(jj,kk))
                T1c1(jj,kk) = 0;
            else
            end
        end
    end

    figure(6);
    
    T1c1min = min(T1c1(T1c1>0)); 
    T1c1max = max(T1c1(:));  
    T1c1_rescaled = imagesc(T1c1, [T1c1min T1c1max]);
    colormap gray

    title({'Figure 6. Final T1c image', 'estimated from FAIR-EPI'})
    name_fig6 = sprintf("%s\\%s_Fig6.png", slicedirname, filename_export); 
    saveas(gcf,name_fig6)

    % Save uncorrected T1 DICOM 
    T1c1_dcm = im2int16(mat2gray(T1c1)); 
    T1c1_dcm_masked = int16(Mask_output).*T1c1_dcm; 
    filename_dcm = sprintf("%s\\%s_T1c.dcm", slicedirname, filename_export); 
    dicomwrite(T1c1_dcm_masked,filename_dcm,q1);

    % Upsample T1c1 image
    T1c1_resized = imresize(T1c1, [pnR rnR]); 

    [optimizer,metric] = imregconfig("multimodal");
    optimizer.InitialRadius = 0.004;

    % Preregister RARE to upsampled T1c  
    ImageR_flipped = fliplr(ImageR);
    ImageR_flipped_resized = imresize(ImageR_flipped, [pn rn]);
    MaskR_flipped = fliplr(MaskR);
    Mask_resized = imresize(Mask, [pnR rnR]); 

    % First, register based on masks / binarized images 
    RARE_tform_step1 = imregtform(MaskR_flipped,Mask_resized,"rigid",optimizer,metric);

    RARE_sameAsInput_step1 = affineOutputView(size(MaskR_flipped),RARE_tform_step1,"BoundsStyle","SameAsInput");
    [MaskR_flipped_step1, MaskR_flipped_step1_ref] = imwarp(MaskR_flipped, RARE_tform_step1, "OutputView", RARE_sameAsInput_step1); 

    RARE_tform_step2 = imregtform(MaskR_flipped_step1,Mask_resized,"similarity",optimizer,metric);

    % Second, apply mask transformation and then register RARE to upsampled T1 from rCBF
    sameAsInput_step1 = affineOutputView(size(ImageR_flipped),RARE_tform_step1,"BoundsStyle","SameAsInput");
    [ImageR_flipped_registered_step1, ImageR_flipped_registered_step1_ref] = imwarp(ImageR_flipped, RARE_tform_step1, "OutputView", sameAsInput_step1); 

    sameAsInput_step2 = affineOutputView(size(ImageR_flipped_registered_step1),RARE_tform_step2,"BoundsStyle","SameAsInput");
    [ImageR_flipped_registered_step2, ImageR_flipped_registered_step2_ref] = imwarp(ImageR_flipped_registered_step1, RARE_tform_step2, "OutputView", sameAsInput_step1); 

    ImageR_flipped_registered_step3 = imregister(ImageR_flipped_registered_step2, T1c1_resized, "similarity",optimizer,metric);
    ImageR_flipped_registered_step4 = imregister(ImageR_flipped_registered_step3, T1c1_resized, "affine",optimizer,metric);

    %% STEP 11: Downsample and register the preregistered RARE reference to T1w image from rCBF fit 
    % First, register based on masks / binarized images 
    ImageR_flipped_registered_resized_step1 = imresize(ImageR_flipped_registered_step4, [pn rn]); 
    %ImageR_flipped_registered_resized_step1 = imresize(ImageR_flipped, [pn rn]); 
    %ImageR_flipped_registered_resized_step1_binarized = double(imbinarize(ImageR_flipped_registered_resized_step1)); 
    ImageR_flipped_registered_resized_step1_binarized = double(ImageR_flipped_registered_resized_step1~=0);

    RARE_resized_tform_step1 = imregtform(ImageR_flipped_registered_resized_step1_binarized,Mask,"rigid",optimizer,metric);

    RARE_resized_sameAsInput_step1 = affineOutputView(size(ImageR_flipped_registered_resized_step1_binarized),RARE_resized_tform_step1,"BoundsStyle","SameAsInput");
    [ImageR_flipped_registered_resized_step1_binarized_step2, ImageR_flipped_registered_resized_step1_binarized_step2_ref] = imwarp(ImageR_flipped_registered_resized_step1_binarized, RARE_resized_tform_step1, "OutputView", RARE_resized_sameAsInput_step1); 

    RARE_resized_tform_step2 = imregtform(ImageR_flipped_registered_resized_step1_binarized_step2,Mask,"similarity",optimizer,metric);

    % Second, apply mask transformation and register the preregistered and downsampled RARE to upsampled T1 from CBF 
    ImageR_flipped_registered_resized_sameAsInput_step1 = affineOutputView(size(ImageR_flipped_registered_resized_step1),RARE_resized_tform_step1,"BoundsStyle","SameAsInput");
    [ImageR_flipped_registered_resized_sameAsInput_step1, ImageR_flipped_registered_resized_sameAsInput_step1_ref] = imwarp(ImageR_flipped_registered_resized_step1, RARE_resized_tform_step1, "OutputView", ImageR_flipped_registered_resized_sameAsInput_step1); 
    ImageR_flipped_registered_resized_sameAsInput_step2 = affineOutputView(size(ImageR_flipped_registered_resized_sameAsInput_step1),RARE_resized_tform_step2,"BoundsStyle","SameAsInput");
    [ImageR_flipped_registered_resized_sameAsInput_step2, ImageR_flipped_registered_resized_sameAsInput_step2_ref] = imwarp(ImageR_flipped_registered_resized_sameAsInput_step1, RARE_resized_tform_step2, "OutputView", ImageR_flipped_registered_resized_sameAsInput_step2); 

    ImageR_flipped_registered_resized_sameAsInput_step3 = imregister(ImageR_flipped_registered_resized_sameAsInput_step2, T1c1, "similarity",optimizer,metric);
    ImageR_flipped_registered_resized = imregister(ImageR_flipped_registered_resized_sameAsInput_step3, T1c1, "affine",optimizer,metric);

    %MaskR_flipped_registered_resized = double(imbinarize(ImageR_flipped_registered_resized)); 
    MaskR_flipped_registered_resized = double(ImageR_flipped_registered_resized~=0); 

    %% STEP 12: Output results  
    figure(7);
    imagesc(ImageR)
    title('Figure 7. Uncorrected RARE anatomical reference')
    name_fig7 = sprintf("%s\\%s_Fig7.png", slicedirname, filename_export); 
    saveas(gcf,name_fig7)

    % Save uncorrected RARE as dicom 
    filename_dcm = sprintf("%s\\%s_RARE_uncorrected.dcm", slicedirname, filename_export); 
    dicomwrite(int16(ImageR),filename_dcm,qR); 

    figure(8);
    imagesc(ImageR_flipped)
    title({'Figure 8. Flipped RARE anatomical reference', 'for registration'})
    name_fig8 = sprintf("%s\\%s_Fig8.png", slicedirname, filename_export); 
    saveas(gcf,name_fig8)

    % Update DCM to reflect new RARE image size 
    qR_new = qR; 
    
    Width_info = dicomfind(qR,"Width");
    Width_info.Value = rn; 
    qR_new = dicomupdate(qR_new,Width_info);
    
    Columns_info = dicomfind(qR,"Columns");
    Columns_info.Value = rn; 
    qR_new = dicomupdate(qR_new,Columns_info);
    
    Height_info = dicomfind(qR,"Height");
    Height_info.Value = pn; 
    qR_new = dicomupdate(qR_new,Height_info);
    
    Rows_info = dicomfind(qR,"Rows");
    Rows_info.Value = pn; 
    qR_new = dicomupdate(qR_new,Rows_info);
    
    % Overlay CBF image on corrected RARE reference without registration 
    imageFused = imfuse(ImageR_flipped_resized, RCBF1); 
    imagesc(imageFused)
    
    figure(9);
    imagesc(imageFused)
    title({'Figure 9. rCBF image', 'overlaid on flipped and resized but unregistered RARE anatomical reference'})
    name_fig9 = sprintf("%s\\%s_Fig9.png", slicedirname, filename_export); 
    saveas(gcf,name_fig9)

    filename_dcm = sprintf("%s\\%s_RARE_flipped_resized.dcm", slicedirname, filename_export); 
    dicomwrite(int16(ImageR_flipped_registered_resized),filename_dcm,qR_new);

    figure(10);
    imageFusedReg = imfuse(ImageR_flipped_registered_resized, RCBF1); 
    imagesc(imageFusedReg)
    title({'Figure 10. rCBF image overlaid on flipped, resized,', ' and registered anatomical reference'})
    name_fig10 = sprintf("%s\\%s_Fig10.png", slicedirname, filename_export); 
    saveas(gcf,name_fig10)
    
    % Save registered RARE DICOM 
    filename_dcm = sprintf("%s\\%s_RARE_flipped_resized_registered.dcm", slicedirname, filename_export); 
    dicomwrite(int16(ImageR_flipped_registered_resized ...
        ),filename_dcm,qR_new);

    % Save FAIR-EPI mask (unprocessed) 
    filename_dcm = sprintf("%s\\FAIR_EPI_mask_origsize_%s.dcm", slicedirname, filename_export); 
    dicomwrite(int16(Mask ...
        ),filename_dcm,q1);

    % Save RARE mask (unprocessed) 
    filename_dcm = sprintf("%s\\RARE_mask_origsize_%s.dcm", slicedirname, filename_export); 
    dicomwrite(int16(MaskR ...
        ),filename_dcm,qR);

    % Save RARE mask (unprocessed) 
    filename_dcm = sprintf("%s\\RARE_mask_resized_%s.dcm", slicedirname, filename_export); 
    dicomwrite(int16(MaskR_flipped_registered_resized ...
        ),filename_dcm,qR_new);

    disp(sprintf('Analysis %d of %d is complete.', dataset_ii, num_datasets)); 

end

close all

end
