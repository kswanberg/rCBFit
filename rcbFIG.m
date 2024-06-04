%function rcbFIG()

%% STEP 1: Define filepaths for masked rCBF statistics calculation 
close all; 
clear all; 

% Create directory for output files 
slicedirroot = pwd; 
datetime_var = strrep(strrep(strrep(string(datetime), ':', '_'), ' ', '_'), '-', '_');
slicedirname = sprintf("%s\\%s_rcbFIGOutput", slicedirroot, datetime_var);

% Prompt user for batch input CSV
uiwait(msgbox('Please select your batch-processing input CSV (Columns: rCBF_Filepath, Mask_Filepath)!'))
[input_CSV_path,~] = uigetfile({'*.csv'}, 'Select your input CSV.', slicedirroot,'MultiSelect','off');
opts = detectImportOptions(input_CSV_path,'ReadVariableNames', false, 'HeaderLines', 1);
%opts = setvartype(opts, 'Mask_Filepath', 'char'); 
opts.Delimiter = ','; 
input_CSV = readtable(input_CSV_path, opts);
num_datasets = height(input_CSV); 

% Check that all folders and files from input CSV actually exist and will not cause problems later  
for dataset_ii = 1:num_datasets

    for inputCSVcol_kk = 1:2

        inputPath = string(input_CSV{dataset_ii, inputCSVcol_kk});
        
        if ~strcmp(inputPath, '1')
            if ~exist(inputPath, 'file')
                error('Initial input check has discovered that %s provided in your input CSV does not exist as written', inputPath)
            end
        end
    end

end 

disp('Initial check of input CSV successful!')

% Initialize vectors for outputs 
rCBF_image_vector = input_CSV{:,1};
mask_vector = input_CSV{:,2};
mask_size_vector = zeros(num_datasets, 1); 
signal_sum_vector = zeros(num_datasets, 1); 
signal_ave_vector = zeros(num_datasets, 1); 
signal_sd_vector = zeros(num_datasets, 1); 
signal_abs_sum_vector = zeros(num_datasets, 1); 
signal_abs_ave_vector = zeros(num_datasets, 1); 
signal_abs_sd_vector = zeros(num_datasets, 1); 

% If export directory does not already exist, make it so 
if not(isfolder(slicedirname))
    mkdir(slicedirname) 
end 

% Now walk through and calculate all values 
for dataset_ii = 1:num_datasets
    
    rCBF_image_filepath = string(input_CSV{dataset_ii, 1});
    rCBF_image = dicomread(rCBF_image_filepath)

    mask_filepath = string(input_CSV{dataset_ii, 2});

    if strcmp(mask_filepath, '1')
    
        mask=rCBF_image~=0

    else
    
        mask = imread(mask_filepath); 
        mask = double(mask~=0); 

    end

    [a,filename_export,c]= fileparts((string(rCBF_image_filepath(1))));
   
    figure(1); 
    fused_image = imfuse(mask, rCBF_image,"blend"); 
    imagesc(fused_image); 

    title({'Figure 1. Input mask overlaid on rCBF image'})
    name_fig1 = sprintf("%s\\%s_Fig1.png", slicedirname, filename_export); 
    saveas(gcf,name_fig1)

    rCBF_image_masked = int16(mask).*rCBF_image; 
    rCBF_image_masked_nonzeros = nonzeros(rCBF_image_masked); 

    mask_size_vector(dataset_ii) = sum(mask(:)); 
    signal_sum_vector(dataset_ii) = sum(rCBF_image_masked_nonzeros(:)); 
    signal_ave_vector(dataset_ii) = mean(rCBF_image_masked_nonzeros(:)); 
    signal_sd_vector(dataset_ii) = std(double(rCBF_image_masked_nonzeros(:))); 
    signal_abs_sum_vector(dataset_ii) = sum(abs(rCBF_image_masked_nonzeros(:))); 
    signal_abs_ave_vector(dataset_ii) = mean(abs(rCBF_image_masked_nonzeros(:))); 
    signal_abs_sd_vector(dataset_ii) = std(double(abs(rCBF_image_masked_nonzeros(:)))); 
end

output_table = table(rCBF_image_vector, mask_vector, mask_size_vector, signal_sum_vector, signal_ave_vector, signal_sd_vector,  signal_abs_sum_vector, signal_abs_ave_vector, signal_abs_sd_vector); 

name_table = sprintf("%s\\%s_rcbFIG_Output_table.csv", slicedirname, datetime_var); 
writetable(output_table, name_table); 

%end
