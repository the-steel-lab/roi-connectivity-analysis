%% roi_connectivity_statistics
if exist('idx_type_label','var')
    fprintf('\n========== Statistics: %s ==========\n', idx_type_label);
end

%%
%% STATISTICAL ANALYSIS: LME on ROI Connectivity
% Fits a 3-way LME with SupParRegion, Task, and Surface as fixed effects
% and Subject as a random intercept. Model formula:
%   Connectivity ~ SupParRegion * Task * Surface + (1|Subject)
%
% Factor definitions
%   SupParRegion : which superior parietal seed was used
%                  'SPPA'  = SupPar         (perception-defined, index 4)
%                  'SPMA'  = SupPar-Mem     (memory-defined,     index 8)
%   Task         : which task defined the target surface ROIs
%                  'Per'   = dynloc (OPA, PPA, MPA)
%                  'Mem'   = imagery (LPMA, VPMA, MPMA)
%   Surface      : topographic subdivision of the target ROIs
%                  'Lateral' = OPA  / LPMA  (roi indices 1 / 5)
%                  'Ventral' = PPA  / VPMA  (roi indices 2 / 6)
%                  'Medial'  = MPA  / MPMA  (roi indices 3 / 7)

% ----- exclude subjects with any NaN in their correlation matrix -----
% A subject flagged earlier (percep_mem_ts_matrix set to NaN) will produce
% an all-NaN slice in mean_corr_mat. Identify and remove them so the LME
% and t-tests only operate on complete cases.
all_roi_cols   = [target_idx_per, target_idx_mem, seed_indices];  % all columns of interest
valid_subj_mask = true(1, length(subjects));
for s = 1:length(subjects)
    slice = mean_corr_mat(all_roi_cols, all_roi_cols, s);
    if any(isnan(slice(:)))
        valid_subj_mask(s) = false;
        fprintf('Excluding subject %s from LME — NaN detected in correlation matrix.\n', subjects{s});
    end
end

subjects_lme = subjects(valid_subj_mask);   % subject list used for stats
mean_corr_lme = mean_corr_mat(:,:,valid_subj_mask);  % filtered corr matrix
nSubj_lme = sum(valid_subj_mask);
fprintf('%d of %d subjects retained for LME analysis.\n', nSubj_lme, length(subjects));

% ----- factor and index definitions ----------------------------------
sup_par_labels = {'SPPA', 'SPMA'};   % seed labels
seed_indices   = [4, 8];             % columns in corrmat / mean_corr_mat

task_labels    = {'Per', 'Mem'};     % task used to define target ROIs
surface_labels = {'Lateral', 'Ventral', 'Medial'};

% Target ROI column indices in mean_corr_mat for each task
%   Per-task : OPA(1)  PPA(2)  MPA(3)
%   Mem-task : LPMA(5) VPMA(6) MPMA(7)
target_idx_per = [1, 2, 3];
target_idx_mem = [5, 6, 7];

% ----- build long-format table ---------------------------------------
nSeeds = numel(seed_indices);
nTasks = numel(task_labels);
nSurfs = numel(surface_labels);
nRows  = nSubj_lme * nSeeds * nTasks * nSurfs;

lme_vars  = {'Subject', 'SupParRegion', 'Task', 'Surface', 'Connectivity'};
data_cell = cell(nRows, numel(lme_vars));

c = 0;
for si = 1:nSeeds
    seed_col   = seed_indices(si);
    seed_label = sup_par_labels{si};

    for ti = 1:nTasks
        if ti == 1
            target_cols = target_idx_per;
        else
            target_cols = target_idx_mem;
        end
        task_label = task_labels{ti};

        for surf_i = 1:nSurfs
            target_col  = target_cols(surf_i);
            surf_label  = surface_labels{surf_i};

            for s = 1:nSubj_lme
                c = c + 1;
                % Fisher z-transformed connectivity averaged across hemispheres
                conn_val = atanh(mean_corr_lme(seed_col, target_col, s));
                data_cell{c, 1} = subjects_lme{s};
                data_cell{c, 2} = seed_label;
                data_cell{c, 3} = task_label;
                data_cell{c, 4} = surf_label;
                data_cell{c, 5} = conn_val;
            end
        end
    end
end

T_conn = cell2table(data_cell, 'VariableNames', lme_vars);

% ----- fit LME -------------------------------------------------------
lme_conn   = fitlme(T_conn, 'Connectivity ~ SupParRegion * Task * Surface + (1|Subject)');
anova_conn = anova(lme_conn);

fprintf('\n=== Connectivity LME (SupParRegion * Task * Surface) ===\n');
disp(anova_conn);

% Store results on the workspace
conn_stats.lme   = lme_conn;
conn_stats.anova = anova_conn;

% ----- follow-up t-tests: SupParRegion (SPPA vs SPMA) ----------------
fprintf('\n--- Follow-up: SupParRegion ---\n');
meanT_region = groupsummary(T_conn, {'Subject','SupParRegion'}, 'mean', 'Connectivity');
meanT_region = sortrows(meanT_region, {'SupParRegion','Subject'});

sppa_vals = meanT_region(strcmp(meanT_region.SupParRegion, 'SPPA'), :);
spma_vals = meanT_region(strcmp(meanT_region.SupParRegion, 'SPMA'), :);

[~, p, ~, stats] = ttest(sppa_vals.mean_Connectivity, spma_vals.mean_Connectivity);
diff_r = sppa_vals.mean_Connectivity - spma_vals.mean_Connectivity;
d_r    = mean(diff_r, 'omitnan') / std(diff_r, 0, 'omitnan');
fprintf('SPPA vs SPMA:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
    stats.df, stats.tstat, p, d_r);

conn_stats.ttest_region = struct( ...
    'region1', 'SPPA', 'region2', 'SPMA', ...
    'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_r);

% ----- follow-up t-tests: Task (Per vs Mem) --------------------------
fprintf('\n--- Follow-up: Task ---\n');
meanT_task = groupsummary(T_conn, {'Subject','Task'}, 'mean', 'Connectivity');
meanT_task = sortrows(meanT_task, {'Task','Subject'});

per_vals = meanT_task(strcmp(meanT_task.Task, 'Per'), :);
mem_vals = meanT_task(strcmp(meanT_task.Task, 'Mem'), :);

[~, p, ~, stats] = ttest(per_vals.mean_Connectivity, mem_vals.mean_Connectivity);
diff_t = per_vals.mean_Connectivity - mem_vals.mean_Connectivity;
d_t    = mean(diff_t, 'omitnan') / std(diff_t, 0, 'omitnan');
fprintf('Per vs Mem:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
    stats.df, stats.tstat, p, d_t);

conn_stats.ttest_task = struct( ...
    'task1', 'Per', 'task2', 'Mem', ...
    'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_t);

% ----- follow-up t-tests: Surface (pairwise) -------------------------
fprintf('\n--- Follow-up: Surface (pairwise) ---\n');
meanT_surf_conn = groupsummary(T_conn, {'Subject','Surface'}, 'mean', 'Connectivity');
meanT_surf_conn = sortrows(meanT_surf_conn, {'Surface','Subject'});

pair_idx = 0;
for i = 1:nSurfs
    for j = i+1:nSurfs
        pair_idx = pair_idx + 1;
        surf_i = surface_labels{i};
        surf_j = surface_labels{j};

        vals_i = meanT_surf_conn(strcmp(meanT_surf_conn.Surface, surf_i), :);
        vals_j = meanT_surf_conn(strcmp(meanT_surf_conn.Surface, surf_j), :);

        [~, p, ~, stats] = ttest(vals_i.mean_Connectivity, vals_j.mean_Connectivity);
        diff_s = vals_i.mean_Connectivity - vals_j.mean_Connectivity;
        d_s    = mean(diff_s, 'omitnan') / std(diff_s, 0, 'omitnan');
        fprintf('%s vs %s:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
            surf_i, surf_j, stats.df, stats.tstat, p, d_s);

        conn_stats.ttest_surface(pair_idx) = struct( ...
            'surf1', surf_i, 'surf2', surf_j, ...
            'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_s);
    end
end

%% STATISTICAL ANALYSIS: LME averaged across Surface (SupParRegion * Task)
% This version first averages connectivity across the three surface subdivisions
% (Lateral, Ventral, Medial) within each SupParRegion × Task cell, then fits a
% simpler 2-way LME:
%   Connectivity ~ SupParRegion * Task + (1|Subject)
%
% This treats Surface as a nuisance dimension and tests whether the two seeds
% differ in their overall sensitivity to how the target set was defined (Per vs Mem).

% ----- average connectivity across surfaces per subject/seed/task ----
nRows_avg  = nSubj_lme * nSeeds * nTasks;
avg_vars   = {'Subject', 'SupParRegion', 'Task', 'Connectivity'};
avg_cell   = cell(nRows_avg, numel(avg_vars));

c = 0;
for si = 1:nSeeds
    seed_col   = seed_indices(si);
    seed_label = sup_par_labels{si};

    for ti = 1:nTasks
        if ti == 1
            target_cols = target_idx_per;
        else
            target_cols = target_idx_mem;
        end
        task_label = task_labels{ti};

        for s = 1:nSubj_lme
            c = c + 1;
            % Fisher z-transform each target ROI then average across surfaces
            conn_vals = atanh(squeeze(mean_corr_lme(seed_col, target_cols, s)));
            avg_conn  = mean(conn_vals, 'omitnan');

            avg_cell{c, 1} = subjects_lme{s};
            avg_cell{c, 2} = seed_label;
            avg_cell{c, 3} = task_label;
            avg_cell{c, 4} = avg_conn;
        end
    end
end

T_conn_avg = cell2table(avg_cell, 'VariableNames', avg_vars);

% ----- fit 2-way LME -------------------------------------------------
lme_conn_avg   = fitlme(T_conn_avg, 'Connectivity ~ SupParRegion * Task + (1|Subject)');
anova_conn_avg = anova(lme_conn_avg);

fprintf('\n=== Connectivity LME — Surface-Averaged (SupParRegion * Task) ===\n');
disp(anova_conn_avg);

conn_stats.lme_avg   = lme_conn_avg;
conn_stats.anova_avg = anova_conn_avg;

% ----- follow-up t-tests: SupParRegion (SPPA vs SPMA) ----------------
fprintf('\n--- Follow-up (surface-averaged): SupParRegion ---\n');
meanT_reg_avg = groupsummary(T_conn_avg, {'Subject','SupParRegion'}, 'mean', 'Connectivity');
meanT_reg_avg = sortrows(meanT_reg_avg, {'SupParRegion','Subject'});

sppa_avg = meanT_reg_avg(strcmp(meanT_reg_avg.SupParRegion, 'SPPA'), :);
spma_avg = meanT_reg_avg(strcmp(meanT_reg_avg.SupParRegion, 'SPMA'), :);

[~, p, ~, stats] = ttest(sppa_avg.mean_Connectivity, spma_avg.mean_Connectivity);
diff_r_avg = sppa_avg.mean_Connectivity - spma_avg.mean_Connectivity;
d_r_avg    = mean(diff_r_avg, 'omitnan') / std(diff_r_avg, 0, 'omitnan');
fprintf('SPPA vs SPMA:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
    stats.df, stats.tstat, p, d_r_avg);

conn_stats.ttest_avg_region = struct( ...
    'region1', 'SPPA', 'region2', 'SPMA', ...
    'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_r_avg);

% ----- follow-up t-tests: Task (Per vs Mem) --------------------------
fprintf('\n--- Follow-up (surface-averaged): Task ---\n');
meanT_task_avg = groupsummary(T_conn_avg, {'Subject','Task'}, 'mean', 'Connectivity');
meanT_task_avg = sortrows(meanT_task_avg, {'Task','Subject'});

per_avg = meanT_task_avg(strcmp(meanT_task_avg.Task, 'Per'), :);
mem_avg = meanT_task_avg(strcmp(meanT_task_avg.Task, 'Mem'), :);

[~, p, ~, stats] = ttest(per_avg.mean_Connectivity, mem_avg.mean_Connectivity);
diff_t_avg = per_avg.mean_Connectivity - mem_avg.mean_Connectivity;
d_t_avg    = mean(diff_t_avg, 'omitnan') / std(diff_t_avg, 0, 'omitnan');
fprintf('Per vs Mem:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
    stats.df, stats.tstat, p, d_t_avg);

conn_stats.ttest_avg_task = struct( ...
    'task1', 'Per', 'task2', 'Mem', ...
    'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_t_avg);

% ----- follow-up t-tests: interaction cells (2x2 simple effects) -----
% Tests each SupParRegion separately for the Per vs Mem contrast, and
% each Task separately for the SPPA vs SPMA contrast.
fprintf('\n--- Follow-up (surface-averaged): Interaction simple effects ---\n');

region_levels = {'SPPA','SPMA'};
task_levels   = {'Per','Mem'};

% Simple effect of Task within each SupParRegion
for ri = 1:2
    reg = region_levels{ri};
    sub_T = T_conn_avg(strcmp(T_conn_avg.SupParRegion, reg), :);
    sub_T = sortrows(sub_T, {'Task','Subject'});

    per_sub = sub_T(strcmp(sub_T.Task, 'Per'), :);
    mem_sub = sub_T(strcmp(sub_T.Task, 'Mem'), :);

    [~, p, ~, stats] = ttest(per_sub.Connectivity, mem_sub.Connectivity);
    diff_si = per_sub.Connectivity - mem_sub.Connectivity;
    d_si    = mean(diff_si, 'omitnan') / std(diff_si, 0, 'omitnan');
    fprintf('%s — Per vs Mem:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
        reg, stats.df, stats.tstat, p, d_si);

    conn_stats.ttest_avg_simple_task(ri) = struct( ...
        'SupParRegion', reg, 'task1', 'Per', 'task2', 'Mem', ...
        'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_si);
    conn_diff(:,ri) = diff_si;
end
%
% a priori t-test on connectivity difference
[~, p, ~, stats] = ttest(conn_diff(:,1),conn_diff(:,2));
    
fprintf('Diff. Per v Mem:SPPA v PPMA — Per vs Mem:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
    stats.df, stats.tstat, p, d_si);
diff_si = conn_diff(:,1)-conn_diff(:,2);
d_si    = mean(diff_si, 'omitnan') / std(diff_si, 0, 'omitnan')
conn_stats.ttest_diff_interaction = struct( ...
        'SupParRegion', 'Diff', 'task1', 'Per', 'task2', 'Mem', ...
        'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_si);
%

% Simple effect of SupParRegion within each Task
for ti = 1:2
    tsk = task_levels{ti};
    sub_T = T_conn_avg(strcmp(T_conn_avg.Task, tsk), :);
    sub_T = sortrows(sub_T, {'SupParRegion','Subject'});

    sppa_sub = sub_T(strcmp(sub_T.SupParRegion, 'SPPA'), :);
    spma_sub = sub_T(strcmp(sub_T.SupParRegion, 'SPMA'), :);

    [~, p, ~, stats] = ttest(sppa_sub.Connectivity, spma_sub.Connectivity);
    diff_si = sppa_sub.Connectivity - spma_sub.Connectivity;
    d_si    = mean(diff_si, 'omitnan') / std(diff_si, 0, 'omitnan');
    fprintf('%s — SPPA vs SPMA:\tt(%d) = %.3f, p = %.3f, d = %.3f\n', ...
        tsk, stats.df, stats.tstat, p, d_si);

    conn_stats.ttest_avg_simple_region(ti) = struct( ...
        'Task', tsk, 'region1', 'SPPA', 'region2', 'SPMA', ...
        'tstat', stats.tstat, 'df', stats.df, 'p', p, 'cohens_d', d_si);
end