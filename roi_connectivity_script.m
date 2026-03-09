%% SETUP AND INITIALIZATION
% This section adds necessary MATLAB toolboxes to the path and loads brain coordinate data

% *** ACTION REQUIRED: Update this path to where AFNI MATLAB tools are installed on your computer ***
% Add AFNI MATLAB tools (used for loading neuroimaging data files)
addpath(genpath('~/Documents/MATLAB/afni_matlab'))
addpath(genpath('~/Documents/MATLAB/utils-2025/'))
% *** ACTION REQUIRED: Update this path to where your ROI files are stored ***
% Define the path where ROI (Region of Interest) files are stored
% If the rois folder is in the same directory as this script, you can leave it as './rois/'
path_to_roi = '../roi_centroid_analysis/output/Dartmouth/com_rois/';

% Define hemisphere labels: 'lh' = left hemisphere, 'rh' = right hemisphere
hemis = {'lh','rh'};

% *** ACTION REQUIRED: Update this path to where your brain surface coordinate files are located ***
% Loop through both hemispheres to load 3D coordinate data
% These coordinates define the positions of vertices on the inflated brain surface
for h = 1:2
    % Load coordinates file for this hemisphere and store in a structure
    % coords.lh will contain left hemisphere coordinates
    % coords.rh will contain right hemisphere coordinates
    % You need the files: lh.inflated.coords.txt and rh.inflated.coords.txt
    coords.(hemis{h}) = load(['../../' hemis{h} '.inflated.coords.txt']);
end
%

%% DEFINE ANALYSIS PARAMETERS
% This section specifies which subjects, ROIs, and ROI versions to analyze

% *** ACTION REQUIRED: Update this list with YOUR subject IDs ***
% List of all subject IDs to include in the analysis (n=12 subjects)
% These IDs should match the subject folder names and ROI file naming convention
subjects = {'scim001','scim002','scim008','scim012','scim014',...
    'scim076','scim0989','scim257','scim297','scim357','scim312','scim995'};

% *** ACTION REQUIRED: Update these ROI numbers if you're analyzing different regions ***
% ROI numbers to analyze (these are specific parcels from a brain atlas)
% 39, 40, 41, 42 correspond to specific brain regions
% These numbers should match the ROI numbers in your ROI filenames
rois = [39 40 41 42];

% *** ACTION REQUIRED: Update these labels if your ROI versions have different names ***
% ROI versions to load:
% one set of ROIs comes from each task
tasks = {'dynloc','imagery'};

COM_size = 300;

%% Initialize data storage structure
% Create a cell array to store all subject data
% Each cell will contain a structure with that subject's ROI and fMRI data

subjectdata = cell(size(subjects));
%%

%% LOAD ROI DATA FOR ALL SUBJECTS
% This nested loop loads ROI vertex indices for each subject, hemisphere, ROI, and version
% The result is a structure containing all ROI definitions for all subjects

for s = 1:length(subjects)  % Loop through each subject
    subject = subjects{s};  % Get current subject ID (e.g., 'scim001')
    %sdata = struct;  % Initialize empty structure to hold this subject's data

    for h = 1:2  % Loop through both hemispheres (left and right)
        hemi = hemis{h};  % Get current hemisphere label ('lh' or 'rh')

        for r = 1:length(rois)  % Loop through each ROI number (39, 40, 41, 42)
            for rv = 1:2  % Loop through each ROI version (per, mim)

                % Get the current ROI version label
                task = tasks{rv};

                % Construct the filename for this specific ROI file
                % Format: subject.version.roinumber.hemisphere.1D.roi
                % Example: '39.lh.1D.roi'
                rfile = [path_to_roi subject '/' task '/' num2str(rois(r)) '.COM.' num2str(COM_size) '.' hemi '.1D.roi'];

                % Display which file is being loaded (helpful for tracking progress)
                disp(['loading rois: ' rfile])

                % Load the ROI file, which contains vertex indices
                % These indices specify which vertices on the brain surface belong to this ROI
                % Store in a nested structure: sdata.per.roi_39.lh.idx (for example)
                try
                    sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx = load(rfile);

                    % remove header if present by checking for ROI size in
                    % first row. this is because of how the COM roi
                    % generation script works
                    if size(sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx,1) == COM_size+1

                        sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx =...
                            sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx(2:COM_size+1,:);
                        disp(size(sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx))
                    end
                    
                catch
                    sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx = NaN;
                    fprintf('subject: %s hemi %s roi %2d task %s does not exist\n', ...
                        subject, hemi, rois(r), task)
                end
            end
        end

    end
    % Store this subject's complete data structure in the main cell array
    subjectdata{s} = sdata;
end

%% LOAD RESTING STATE fMRI DATA
% This section loads preprocessed resting state (rs) functional MRI data for each subject
% Resting state data shows brain activity when the subject is not performing a task

for s = 1:length(subjects)  % Loop through each subject

    subject = subjects{s};  % Get current subject ID
    sdata = subjectdata{s};  % Get this subject's data structure
    disp(subject)  % Display subject ID to track progress

    for h = 1:2  % Loop through both hemispheres
        hemi = hemis{h};  % Get current hemisphere label

        % Find the directory containing this subject's preprocessed data
        % The '*manual' wildcard matches a directory with 'manual' in the name
        % dir() returns information about the matching directory
        path_to_rsdata = dir(['~/Dropbox/Projects/social-imagery/data/' subject '/ppTedAna/*manual']);

        % Reconstruct the full path to the data directory
        % Combines the folder path and the directory name
        path_to_rsdata = [path_to_rsdata.folder '/' path_to_rsdata.name '/'];

        % Construct the filename for the resting state data
        % This data has been:
        %   - denoised (noise removed)
        %   - GSR (global signal regression applied)
        %   - z-scored (normalized)
        %   - smoothed (spatially blurred to reduce noise)
        rsdata_file = [path_to_rsdata 'desc-denoised_bold.GSR.zscore.smooth.' hemi '.1D.dset'];

        % Load the resting state data using AFNI's BrikLoad function
        % This creates a matrix where rows are vertices and columns are time points
        sdata.rsdata.(hemi) = BrikLoad(rsdata_file);
    end

    % Update the main data structure with the resting state data
    subjectdata{s} = sdata;
end


%% EXTRACT ROI TIME SERIES FROM CENTER OF MASS ROIs
% This section extracts the average time series from each COM ROI
% A time series shows how activity in that brain region changes over time

% Initialize a 4D matrix to store all time series data
% Dimensions: 300 time points × 8 ROIs × 2 hemispheres × 12 subjects
% We use 8 ROIs because we have 4 ROI numbers × 2 versions (per and mim, excluding comb)
percep_mem_ts_matrix = nan(300,8,2,length(subjects));

for s = 1:length(subjects)  % Loop through each subject
    sdata = subjectdata{s};  % Get this subject's data structure
    subject = subjects{s};  % Get subject ID
    disp(subject)  % Display subject ID to track progress

    for h = 1:2  % Loop through both hemispheres
        hemi = hemis{h};  % Get current hemisphere label
        rcount = 0;  % Initialize counter for tracking which ROI column we're filling

        for rv = 1:2  % Loop through ROI versions (per, mim, comb)
            for r = 1:length(rois)  % Loop through each ROI number
                task = tasks{rv};  % Get ROI version label

                % Get the COM vertex indices for this ROI
                current_roi_indices = sdata.(task).(['roi_' num2str(rois(r))]).(hemi).idx;
                
                % Extract the time series for this ROI
                % This averages the activity across all vertices in the ROI at each time point
                % The first output (~) is discarded, the second output is the time series
                [~, roi_ts] = extract_roi_timeseries(sdata.rsdata.(hemi),current_roi_indices);

                % Store the time series in the data structure
                sdata.(task).(['roi_' num2str(rois(r))]).(hemi).com_rs_ts = roi_ts;


                    rcount = rcount+1;  % Increment ROI counter

                    % Store this ROI's time series in the master matrix
                    % This creates a matrix with all subjects' data for statistical analysis
                    percep_mem_ts_matrix(:,rcount,h,s) = roi_ts;

                    % Create a label for this ROI (for plotting later)
                    roi_labels{rcount} = [ task ' roi ' num2str(rois(r)) ];
                
            end
        end
    end
        subjectdata{s} = sdata;

end
%% filter corrmat to remove any subjects with any nans in their matrix
for s = 1:12
    for h = 1:2
        if any(isnan(percep_mem_ts_matrix(:,:,h,s)),'all')
            disp(['NaN detected for : ' subjects{s} ' hemisphere ' hemis{h}]);
            percep_mem_ts_matrix(:,:,h,s) = nan(size(percep_mem_ts_matrix(:,:,h,s)));
        end
    end
end
%


%% CALCULATE CORRELATIONS AMONG ROI TIME SERIES
% This computes functional connectivity between all pairs of ROIs
% Correlation measures how similarly two brain regions' activity patterns change over time
% High correlation = strong functional connectivity = regions work together

for s = 1:12  % Loop through all 12 subjects
    for h = 1:2  % Loop through both hemispheres

        % Calculate pairwise correlations between all ROIs' time series
        % Input: 300 timepoints × 8 ROIs matrix
        % Output: 8×8 correlation matrix
        % Each cell (i,j) shows how strongly ROI i correlates with ROI j
        % The diagonal will be 1.0 (each ROI perfectly correlates with itself)
        corrmat(:,:,h,s) = corr(percep_mem_ts_matrix(:,:,h,s));
    end
end

%% PLOT CONNECTIVITY DATA FROM ALL REGIONS
% This section creates a figure showing connectivity patterns between ROIs
% The analysis focuses on two "seed" regions and their connectivity to other regions

% Define human-readable names for the 8 ROIs
% These correspond to the 8 columns in percep_mem_ts_matrix
roi_names = {'OPA','PPA','MPA','SupPar','LPMA','VPMA','MPMA','SupPar-Mem'};

% Define which ROIs are the "seeds" (main regions of interest)
% Index 4 = SupPar (perception), Index 8 = SupPar-Mem (memory)
bars_of_interest = [4 8];

% Create a logical index to exclude the seed ROIs themselves from the plot
% We want to show connectivity TO other regions, not self-correlation
baridx = bar_order;
baridx = ~ismember(baridx,bars_of_interest);

% Define the order in which to display ROIs in the plots
% This reorders the ROIs for better visualization
bar_order = [1 5 2 6 3 7 4 8];

% Average the correlation matrix across hemispheres
% squeeze() removes the singleton dimension (hemisphere) after averaging
% Result: 8×8×12 matrix (ROIs × ROIs × subjects)
mean_corr_mat = squeeze(mean(corrmat,3,'omitnan'));

% Define indices for grouping perception vs memory ROIs
% Row 1: perception ROIs (indices 1,2,3)
% Row 2: memory ROIs (indices 5,6,7)
sub_inds = [1 2 3; 5 6 7];

% Create a new figure window
figure()
YLIM = [-0.5 1.5];
% Loop through the two seed ROIs
for r = 1:2

    % FIRST SET OF SUBPLOTS: Bar plots showing connectivity to all other regions
    % Position subplots in locations specified by sub_inds
    subplot(2,4,sub_inds(r,:),'linewidth',2)
    hold on  % Allow multiple plot elements to be drawn on same axes

    % Extract correlation data for this seed ROI
    % atanh() performs Fisher's z-transformation to normalize correlations
    % This transformation makes correlations more suitable for statistical analysis
    bardata = atanh(mean_corr_mat(bars_of_interest(r),:,:));

    % Plot bar graph showing mean connectivity across subjects
    % Dimension 3 is subjects, so mean(...,3) averages across all 12 subjects
    bar(mean(bardata(:,bar_order(baridx),:),3,'omitnan'),'barwidth',0.75);

    % Overlay individual subject data as connected points
    % This shows individual variability around the group mean
    plot(squeeze(bardata(:,bar_order(baridx),:)),'-k','linewidth',1)
    plot(squeeze(bardata(:,bar_order(baridx),:)),'ok','markerfacecolor','k','MarkerEdgeColor','none')

    % Set x-axis tick positions (one for each non-seed ROI)
    xticks([1:6])
    xlim([0.5 6.5])

    % Label x-axis with ROI names (excluding the seed ROIs)
    xticklabels(roi_names(bar_order(baridx)))

    % Set y-axis limits to show the range of correlations
    ylim(YLIM)
    ylabel('Correlation (z)')  % Label y-axis (z = Fisher-transformed)
    yticks(YLIM)  % Set y-axis tick positions


    % SECOND SET OF SUBPLOTS: Compare perception vs memory connectivity
    % Position subplot at the seed ROI's index
    subplot(2,4,bars_of_interest(r),'linewidth',2)
    hold on

    % Calculate mean connectivity to perception ROIs (indices 1,2,3)
    % First dimension of cond_means = condition (1=perception, 2=memory)
    cond_means(1,:,r) = squeeze(mean(bardata(:,sub_inds(1,:),:),2,'omitnan'));

    % Calculate mean connectivity to memory ROIs (indices 5,6,7)
    cond_means(2,:,r) = squeeze(mean(bardata(:,sub_inds(2,:),:),2,'omitnan'));

    % Plot bar graph comparing mean perception vs memory connectivity
    bar(mean(cond_means(:,:,r),2,'omitnan'),'barwidth',0.75);

    % Overlay individual subject data as connected points
    % Each line connects one subject's perception and memory connectivity values
    plot(cond_means(:,:,r),'k-','linewidth',1)
    plot(cond_means(:,:,r),'ko','markerfacecolor','k','MarkerEdgeColor','none')

    % Set up x-axis: two bars for perception and memory
    xticks([1 2])
    xticklabels({'Perc','Mem'})
    xlim([.5 2.5])  % Add padding around the two bars

    % Match y-axis formatting to the other subplots
    ylim(YLIM)
    ylabel('Correlation (z)')
    yticks(YLIM)
end

%% run the stats
roi_connectivity_statistics
%%
% for each subject and each hemisphere, correlate each ROI time series with the rest of the brain and save the
% output
path_for_corrmaps = './corrmaps/';
percentile_thresh = 95;
for si = 1:numel(subjects)
    sdata = subjectdata{si};
    subject = subjects{si};

    for hi = 1:2
        hemi = hemis{hi};
        braints = sdata.rsdata.(hemi);

        % correlate roi time series with the rest of the brain
        for ri = 1:4
            roi = rois(ri);
            
            for ti = 1:2
                task = tasks{ti};

                roits = sdata.(task).(['roi_' num2str(roi)]).(hemi).com_rs_ts;
            
                nodeidx = braints(:,1);
            % Perform correlation between ROI time series and brain data
                corr_results = corr(roits', braints(:,2:end)', 'Rows', 'complete');
                % replace nan with 0
                corr_results(isnan(corr_results)) = 0;  % Replace NaNs with 0
                disp([subject '.' task '.roi_' num2str(roi) '.' hemi ' correlated'])

                % create thresholded version of corr_result that keeps only
                % values greater than percentile treshold
                % Apply threshold to correlation results

                thresholded_corr = corr_results > prctile(corr_results, percentile_thresh);
            % save as brain data out as a matrix
                writematrix([nodeidx corr_results'],[path_for_corrmaps '/' task '/' subject '.' task '.roi_' num2str(roi) '.' hemi '.1D.dset'],'FileType','text')
                writematrix([nodeidx thresholded_corr'],[path_for_corrmaps '/' task '/' subject '.' task '.roi_' num2str(roi) '.' num2str(percentile_thresh) 'percentile' hemi '.1D.dset'],'FileType','text')
                
            end
        end
    end
end






