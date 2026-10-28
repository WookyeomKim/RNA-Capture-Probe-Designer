%% Rendering of the genome-wide accessibility landscape
% Step 2 of 2. Produces Supplementary Figure S4C from the data computed by
% P_profile_T_Sweep.m.
%
% WHAT IT DOES
%   Loads the precomputed accessibility landscape and renders four figures:
%     Figure 1  P_unpaired along the HIV-1 genome at 25 and 65 C, 70-nt
%               windows. This is the published Supplementary Figure S4C.
%     Figure 2  Genome-wide mean P_unpaired against temperature, 70-nt only.
%     Figure 3  The same trend for 19-nt and 70-nt windows.
%     Figure 4  The same trend for all four probe lengths.
%   Section 3 then reads off P_unpaired and its genome-wide percentile rank at
%   a set of specified genome positions, and writes them to an Excel file.
%
% REQUIREMENTS
%   MATLAB R2020b or later (uses exportgraphics).
%   HIV_Thermal_Data_Saved2.mat, produced by P_profile_T_Sweep.m. That script
%   must be run first; this one does no folding of its own and will fail
%   immediately if the .mat file is absent.
%   Section 3 additionally writes .xlsx, which needs no extra toolbox.
%
% HOW TO RUN
%   Put HIV_Thermal_Data_Saved2.mat in the working directory, open this file
%   in MATLAB and press Run, or call it by name from the command window.
%   No input arguments are required. Rendering takes a few seconds.
%   Outputs: four PNGs at 300 dpi, named in the Filename1-4 variables below,
%   and one .xlsx table from section 3.
%
% NOTE ON UNITS
%   P_unpaired is a probability on a 0 to 1 scale, not a percentage. The
%   percentile ranks reported in section 3 are percentages.

clc;
clear;
close all;

%% ===================== 1. LOAD DATA =====================
% Loads the precomputed landscape so the long scan does not have to be rerun.
load('HIV_Thermal_Data_Saved2.mat');
fprintf('Data loaded. Rendering figures...\n');

Filename1 = 'P_profile_Across_HIV1.png';
Filename2 = 'P_profile_Across_HIV2.png';
Filename3 = 'P_profile_Across_HIV3.png';
Filename4 = 'P_profile_Across_HIV4.png';

%% ===================== 2. PLOTTING =====================
trend_colors = lines(4);
idx_70 = find(probeLength_range == 70, 1);
idx_19 = find(probeLength_range == 19, 1);

% -----------------------------------------------------------------
% [Figure 1] Genome-wide landscape (70 nt, two representative temperatures)
% This is the published Supplementary Figure S4C.
% -----------------------------------------------------------------
colors_show = jet(num_T_show);
fig1 = figure('Name', 'Fig 1: Landscape (70nt)');
box on; grid on; hold on;
P_land_70 = P_landscape_all{idx_70};

for t_idx = 1:num_T_show
    real_idx = find(T_sweep == T_show(t_idx));
    if ~isempty(real_idx)
        % Plotted as a continuous line, without point markers, for legibility
        plot(scan_indices, P_land_70(real_idx, :), '-', 'Color', colors_show(t_idx, :), ...
            'LineWidth', 2.5, 'DisplayName', sprintf('%d \\circC', T_show(t_idx)));
    end
end
set(gca, 'FontSize', 16);
xlabel('Position on HIV Sequence (nt)');
ylabel('Unpaired Prob');
legend('Location', 'northwest', 'FontSize', 14);
xlim([1 actual_limit]);
exportgraphics(gcf, Filename1, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 2] Global temperature trend (70 nt only)
% -----------------------------------------------------------------
fig2 = figure('Name', 'Fig 2: Global Trend (70nt)');
box on; grid on; hold on;
global_mean = mean(P_land_70, 2, 'omitnan')';
global_std  = std(P_land_70, 0, 2, 'omitnan')';
upper_bound = min(1, global_mean + global_std);
lower_bound = max(0, global_mean - global_std);

X_fill = [T_sweep, fliplr(T_sweep)];
Y_fill = [upper_bound, fliplr(lower_bound)];
fill(X_fill, Y_fill, trend_colors(3,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(T_sweep, global_mean, '-', 'Color', trend_colors(3,:), 'LineWidth', 3.5, ...
    'DisplayName', '70nt (Mean \pm 1 SD)');
set(gca, 'FontSize', 16);
xlabel('Temperature (\circC)');
ylabel('Global Unpaired Prob');
legend('Location', 'northwest', 'FontSize', 16);
xticks(T_min:10:T_max);
xlim([25 85]);
exportgraphics(gcf, Filename2, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 3] Global temperature trend, 19 nt vs 70 nt
% -----------------------------------------------------------------
fig3 = figure('Name', 'Fig 3: Trend (19nt vs 70nt)');
box on; grid on; hold on;
compare_idx = [idx_19, idx_70];
comp_colors = [trend_colors(1,:); trend_colors(3,:)];

% Shaded SD bands first
for k = 1:2
    c_idx = compare_idx(k);
    P_mat = P_landscape_all{c_idx};

    g_mean = mean(P_mat, 2, 'omitnan')';
    g_std  = std(P_mat, 0, 2, 'omitnan')';

    u_b = min(1, g_mean + g_std);
    l_b = max(0, g_mean - g_std);

    fill([T_sweep, fliplr(T_sweep)], [u_b, fliplr(l_b)], comp_colors(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
% Mean lines on top
for k = 1:2
    c_idx = compare_idx(k);
    P_mat = P_landscape_all{c_idx};
    g_mean = mean(P_mat, 2, 'omitnan')';
    plot(T_sweep, g_mean, '-', 'Color', comp_colors(k,:), 'LineWidth', 3.5, ...
        'DisplayName', sprintf('%dnt (Mean \\pm 1 SD)', probeLength_range(c_idx)));
end
set(gca, 'FontSize', 16);
xlabel('Temperature (\circC)');
ylabel('Global Unpaired Prob');
legend('Location', 'northwest', 'FontSize', 14);
xticks(T_min:10:T_max);
xlim([25 85]);
exportgraphics(gcf, Filename3, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 4] Global temperature trend, all four probe lengths
% -----------------------------------------------------------------
fig4 = figure('Name', 'Fig 4: Trend (All Lengths)');
box on; grid on; hold on;

% 1) Shaded SD bands first
for k = 1:num_cases
    P_mat = P_landscape_all{k};
    g_mean = mean(P_mat, 2, 'omitnan')';
    g_std  = std(P_mat, 0, 2, 'omitnan')';

    u_b = min(1, g_mean + g_std);
    l_b = max(0, g_mean - g_std);

    fill([T_sweep, fliplr(T_sweep)], [u_b, fliplr(l_b)], trend_colors(k,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
% 2) Mean lines drawn on top
for k = 1:num_cases
    P_mat = P_landscape_all{k};
    g_mean = mean(P_mat, 2, 'omitnan')';
    plot(T_sweep, g_mean, '-', 'Color', trend_colors(k,:), 'LineWidth', 3.5, ...
        'DisplayName', sprintf('%dnt (Mean \\pm 1 SD)', probeLength_range(k)));
end
set(gca, 'FontSize', 16);
xlabel('Temperature (\circC)');
ylabel('Global Unpaired Prob');
legend('Location', 'northwest', 'FontSize', 14);
xticks(T_min:10:T_max);
xlim([25 85]);
exportgraphics(gcf, Filename4, 'Resolution', 300);

fprintf('All images (Fig 1-4) exported.\n');

%% ===================== 3. SPECIFIC PROBE REGIONS EVALUATION =====================
% Reads off P_unpaired and its genome-wide percentile rank at a set of
% specified genome positions, using the 70-nt landscape at a single
% temperature.
target_T = 65; % evaluation temperature (C)

% Genome positions to evaluate (column vector), given as the start position
% of the 70-nt capture window on the consensus sequence
target_positions = [11031; 6376; 2196; 10971; 3768; 5676; 13005];

% Index of the evaluation temperature within T_sweep
idx_65 = find(T_sweep == target_T);

if isempty(idx_65)
    fprintf('\nWarning: %d C is not present in T_sweep.\n', target_T);
else
    fprintf('\nEvaluating P_unpaired and percentile rank at %d positions (%d C)...\n', length(target_positions), target_T);

    % Full 70-nt probability profile at the evaluation temperature
    P_profile = P_land_70(idx_65, :);

    % Drop NaNs so that interpolation and percentile ranking use real values only
    valid_idx = ~isnan(P_profile);
    valid_P = P_profile(valid_idx);
    valid_scan = scan_indices(valid_idx);

    % The landscape is sampled every 10 nt, so positions in between are
    % obtained by linear interpolation
    target_P_unpaired = interp1(valid_scan, valid_P, target_positions, 'linear');

    % Percentile rank: the fraction of all scanned windows whose P_unpaired is
    % at or below this one. Higher values indicate greater predicted
    % accessibility.
    target_percentiles = zeros(length(target_positions), 1);
    for i = 1:length(target_positions)
        target_percentiles(i) = (sum(valid_P <= target_P_unpaired(i)) / length(valid_P)) * 100;
    end

    ExportTable = table(target_positions, target_P_unpaired, target_percentiles, ...
        'VariableNames', {'Target_Position_nt', 'P_unpaired_Value', 'Percentile_Rank_Pct'});

    disp('--------------------------------------------------------------');
    disp(ExportTable);
    disp('--------------------------------------------------------------');

    excel_filename = sprintf('Algorithm_Targets_Evaluation_70nt_%dC.xlsx', target_T);
    writetable(ExportTable, excel_filename);

    fprintf('Evaluation complete.\n');
    fprintf('Saved: %s\n\n', excel_filename);
end
