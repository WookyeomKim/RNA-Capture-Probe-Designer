% Thermodynamic profiling of the 266-nt synthetic RNA target.
% Compares target accessibility (ViennaRNA/RNAplfold) with probe-target
% hybridization fraction (Sugimoto 1995 RNA/DNA nearest-neighbor model)
% across temperature, for capture regions of 19, 40, 70 and 120 nt.
% Generates Supplementary Figure S4A.
%
% REQUIREMENTS
%   MATLAB R2020b or later (uses readlines and exportgraphics).
%   ViennaRNA 2.6.4 or later, installed locally. RNAplfold is called as an
%   external executable, so the ViennaRNA binaries must be present on disk.
%   Download: https://www.tbi.univie.ac.at/RNA/
%
% SETUP
%   By default the script expects RNAplfold at:
%       ~/viennarna/bin/RNAplfold
%   If ViennaRNA is installed elsewhere, edit the rnaPlfold path inside the
%   function DGopen_for_probe_on_target_RNA near the bottom of this file.
%   On Windows, point it to RNAplfold.exe.
%
% HOW TO RUN
%   Open this file in MATLAB and press Run, or call it by name from the
%   command window. No input arguments are required; all sequences and
%   settings are defined in sections 1 and 2 below. Runtime is a few minutes.
%   The script writes three PNG files (300 dpi) to the working directory,
%   named in the Filename1-3 variables below.

clc;
clear;
close all;

Filename1='Tm_Fraction_vs_Pro_10nt_avg_DnaRnaPar.png';
Filename2='Tm_Fraction_vs_Pro2_10nt_avg_DnaRnaPar.png';
Filename3='Tm_Fraction_vs_Pro3_10nt_avg_DnaRnaPar.png';

%% ===================== 1. INPUT SEQUENCES =====================
% Full 266-nt in vitro transcript (T7 +1 G followed by the gBlock template
% downstream of the promoter). Used as the folding context for RNAplfold.
target = 'GCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGCGGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCT';

% Capture probe binding sites within the target (RNA sense)
CP_19nt  = 'GCCGATAGACAAGGAACTG';
CP_40nt  = 'AGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAA';
CP_70nt  = 'CTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCT';
CP_120nt = 'GGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCT';

cap_targets = {CP_19nt, CP_40nt, CP_70nt, CP_120nt};
cap_names   = {'19nt', '40nt', '70nt', '120nt'};
num_cases   = length(cap_targets);
x_lens      = [19, 40, 70, 120]; % x-axis values for Figure 3

%% ===================== 2. ANALYSIS SETTINGS =====================
targetExtension = 300;   % flanking sequence retained around the binding site

% Experimental strand concentrations
Acon = 1.7e-17;  % RNA target
Bcon = 4e-8;     % DNA capture probe

unpair_window   = 10;    % length of the contiguous unpaired stretch scored

% Statistic used to summarize P_unpaired across the binding site
stat_mode       = 'avg';  % 'avg', 'max', or 'topX'
top_percent     = 0.10;   % used only when stat_mode is 'topX'
salt_M          = 1.0;    % monovalent cation concentration (M)

T_min = 25;
T_max = 100;
dT    = 5;
T_list = T_min:dT:T_max;
T_fine = T_min:0.1:T_max;

%% ===================== 3. DATA COLLECTION LOOP =====================
% Storage for all conditions (used to share y-axis limits and to build
% Figures 2 and 3)
F_fine_all    = cell(num_cases, 1);
P_fine_all    = cell(num_cases, 1);
F_nupack_all  = cell(num_cases, 1);
P_vienna_all  = cell(num_cases, 1);
Tm_open_all   = zeros(num_cases, 1);
P_at_Tm_all   = zeros(num_cases, 1);
Tm_nupack_all = zeros(num_cases, 1); % variable name kept; values are Sugimoto 1995
max_P_global = 0; % global maximum, used to equalize the probability axis

fprintf('Starting calculations for %d conditions...\n', num_cases);
fprintf('Statistic: %s (mean of top %d%%)\n', stat_mode, round(top_percent * 100));
fprintf('Concentrations: target = %g M, probe = %g M\n\n', Acon, Bcon);

for i = 1:num_cases
    cap_tgt = cap_targets{i};
    probe = dna_or_rna_to_rna_revcomp(cap_tgt);

    fprintf('[%d/%d] Processing capture region %s...\n', i, num_cases, cap_names{i});

    % 1. RNA-DNA hybridization fraction from Sugimoto 1995 NN parameters
    [Tm_nupack, T_nupack_list, F_nupack] = Tm_from_Sugimoto1995_RNADNA(cap_tgt, Acon, Bcon, T_min, T_max, dT);

    % 2. ViennaRNA: probability that a contiguous unpair_window-nt stretch is unpaired
    P_vienna = nan(size(T_list));
    for j = 1:numel(T_list)
        tempC = T_list(j);
        out_v = DGopen_for_probe_on_target_RNA(target, probe, tempC, salt_M, targetExtension, unpair_window, stat_mode, top_percent);
        P_vienna(j) = out_v.Pacc;
    end

    % Interpolation for smooth curves
    P_fine = pchip(T_list, P_vienna, T_fine);
    F_fine = pchip(T_list, F_nupack, T_fine);

    % Temperature of maximum opening rate (steepest rise in P_unpaired)
    dP_dT = gradient(P_fine, T_fine);
    [~, max_idx] = max(dP_dT);

    % Store results
    F_fine_all{i}   = F_fine;
    P_fine_all{i}   = P_fine;
    F_nupack_all{i} = F_nupack;
    P_vienna_all{i} = P_vienna;
    Tm_open_all(i)  = T_fine(max_idx);
    P_at_Tm_all(i)  = P_fine(max_idx);
    Tm_nupack_all(i)= Tm_nupack;

    % Update global maximum
    max_P_global = max(max_P_global, max(P_fine));
end

fprintf('\n=========================================================\n');
fprintf('        [Thermodynamic summary: Figures 2 and 3]\n');
fprintf('---------------------------------------------------------\n');
for i = 1:num_cases
    fprintf(' %5s | Hybridization Tm = %5.1f deg C  |  T_open = %5.1f deg C\n', ...
        cap_names{i}, Tm_nupack_all(i), Tm_open_all(i));
end
fprintf('=========================================================\n\n');
fprintf('Calculations complete. Generating figures...\n\n');

%% ===================== 4. PLOTTING =====================
% -----------------------------------------------------------------
% [Figure 1] 2x2 subplots, shared y-axis limits across all four panels
% -----------------------------------------------------------------
fig1 = figure('Name', 'Hybridization vs ViennaRNA (Subplots)', 'Position', [50, 50, 1200, 900]);
for i = 1:num_cases
    subplot(2, 2, i);
    box on; grid on; hold on;

    % Left axis: hybridization fraction
    yyaxis left;
    plot(T_fine, F_fine_all{i}, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5, 'DisplayName', 'Hybrid Fraction');
    plot(T_list, F_nupack_all{i}, 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    ylabel('Hybridization Fraction');
    ylim([0 1.05]);
    ax = gca; ax.YColor = [0 0.4470 0.7410];

    % Right axis: unpaired probability
    yyaxis right;
    plot(T_fine, P_fine_all{i}, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5, 'DisplayName', 'Unpair Prob');
    plot(T_list, P_vienna_all{i}, 'o', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');

    ylabel('Unpaired Probability');
    ylim([0, max(max_P_global * 1.1, 0.01)]);
    ax = gca; ax.YColor = [0.8500 0.3250 0.0980];

    set(gca, 'FontSize', 16);
    title(sprintf('Target Region: %s', cap_names{i}), 'FontSize', 18);
    xlabel('Temperature (\circC)');
    if i == 1, legend('Location', 'best'); end
end
exportgraphics(gcf, Filename1, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 2] Combined plot, 19-nt vs 40-nt only
% -----------------------------------------------------------------
fig2 = figure('Name', 'Combined Plot (19nt vs 40nt)');
box on; grid on; hold on;
colors = lines(num_cases);
h_lines = gobjects(2, 1); % handles for the 19-nt and 40-nt traces

% Left axis: hybridization fraction
yyaxis left;
ylabel('Hybridization Fraction');
ylim([0 1.05]);
ax = gca; ax.YColor = 'k';
for i = 1:2 % plot only case 1 (19 nt) and case 2 (40 nt)
    h_lines(i) = plot(T_fine, F_fine_all{i}, '-', 'Color', colors(i,:), 'LineWidth', 3, 'DisplayName', cap_names{i});
end

% Right axis: unpaired probability
yyaxis right;
ylabel('Unpaired Probability');
ylim([0, max(max_P_global * 1.1, 0.01)]);
ax = gca; ax.YColor = 'k';
for i = 1:2 % plot only case 1 (19 nt) and case 2 (40 nt)
    plot(T_fine, P_fine_all{i}, '--', 'Color', colors(i,:), 'LineWidth', 3, 'HandleVisibility', 'off');
end

% Dummy lines used only to build the legend
h_dummy_frac = plot(nan, nan, 'k-', 'LineWidth', 3, 'DisplayName', 'Hyb. Frac.');
h_dummy_prob = plot(nan, nan, 'k--', 'LineWidth', 3, 'DisplayName', 'Unpair Prob.');

set(gca, 'FontSize', 16);
xlabel('Temperature (\circC)');
%title('Comparison: 19nt vs 40nt', 'FontSize', 18);
h_leg = legend([h_lines(1:2); h_dummy_frac; h_dummy_prob]);
xlim([25 100]);
set(h_leg, 'Position', [0.2, 0.64, 0.15, 0.2]);
exportgraphics(gcf, Filename2, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 3] Delta T (thermodynamic gap) trend across all four lengths
% -----------------------------------------------------------------
fig3 = figure('Name', 'Delta T Trend', 'Position', [150, 150, 900, 600]);
box on; grid on; hold on;
x_fill = [x_lens, fliplr(x_lens)];
y_fill = [Tm_nupack_all', fliplr(Tm_open_all')];

fill(x_fill, y_fill, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', '\Delta T Gap');
plot(x_lens, Tm_nupack_all, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'DisplayName', 'Sugimoto 1995 T_m (Hybridization)');
plot(x_lens, Tm_open_all, '-s', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'DisplayName', 'ViennaRNA T_{open} (Nucleation)');

set(gca, 'FontSize', 16);
xticks(x_lens);
xticklabels({'19nt', '40nt', '70nt', '120nt'});
xlabel('Capture Region Length');
ylabel('Temperature (\circC)');
title('Thermodynamic Gap (\Delta T) vs. Target Length', 'FontSize', 18);
legend('Location', 'best');
exportgraphics(gcf, Filename3, 'Resolution', 300);

%% ===================== FUNCTIONS =====================

function [Tm_C, T_list, frac] = Tm_from_Sugimoto1995_RNADNA(rna_seq, A0, B0, Tmin_C, Tmax_C, dT)
    % Bimolecular hybridization fraction between an RNA target region and its
    % complementary DNA probe, from the RNA/DNA nearest-neighbor parameters of
    % Sugimoto et al. (1995) at 1 M Na+.
    rna_seq = upper(strrep(rna_seq, '-', ''));
    rna_seq = strrep(rna_seq, 'T', 'U');
    N = length(rna_seq);

    % Sugimoto 1995 nearest-neighbor parameters (1 M NaCl)
    NN_dH = containers.Map(...
        {'AA','AC','AG','AU','CA','CC','CG','CU','GA','GC','GG','GU','UA','UC','UG','UU'}, ...
        [-7.8, -5.9, -9.1, -8.3, -9.0, -9.3, -16.3, -7.0, -5.5, -8.0, -12.8, -7.8, -7.8, -8.6, -10.4, -11.5] * 1000);
    NN_dS = containers.Map(...
        {'AA','AC','AG','AU','CA','CC','CG','CU','GA','GC','GG','GU','UA','UC','UG','UU'}, ...
        [-21.9, -12.3, -23.5, -23.9, -26.1, -23.2, -47.1, -19.7, -13.5, -17.1, -31.9, -21.6, -23.2, -22.9, -28.4, -36.4]);

    % Helix initiation terms
    dH_total = 1.9 * 1000;
    dS_total = -3.9;

    for idx = 1:(N-1)
        dinuc = rna_seq(idx:idx+1);
        if isKey(NN_dH, dinuc)
            dH_total = dH_total + NN_dH(dinuc);
            dS_total = dS_total + NN_dS(dinuc);
        else
            error('Invalid dinucleotide found in RNA sequence: %s', dinuc);
        end
    end

    T_list = Tmin_C:dT:Tmax_C;
    frac = nan(size(T_list));
    R = 1.987;

    limiting_con = min(A0, B0);

    for k = 1:length(T_list)
        T_K = T_list(k) + 273.15;
        dG = dH_total - T_K * dS_total;
        K_T = exp(-dG / (R * T_K));

        % 1:1 mass balance, solved for the equilibrium duplex concentration
        S = A0 + B0 + (1 / K_T);
        P = A0 * B0;

        AB_real = (2 * P) / (S + sqrt(S^2 - 4 * P));
        frac(k) = AB_real / limiting_con;
    end

    % Tm taken as the temperature at which the bound fraction crosses 0.5
    idx_half = find(frac >= 0.5, 1, 'last');
    if isempty(idx_half) || idx_half == length(T_list)
        [~, j] = min(abs(frac - 0.5));
        Tm_C = T_list(j);
    else
        T1 = T_list(idx_half); T2 = T_list(idx_half+1);
        f1 = frac(idx_half);   f2 = frac(idx_half+1);
        Tm_C = T1 + (0.5 - f1) * (T2 - T1) / (f2 - f1);
    end
end

function out = DGopen_for_probe_on_target_RNA(target_seq, probe_seq, tempC, salt_M, targetExt, winSize, stat_mode, top_percent)
    % Runs RNAplfold on the target with flanking context and returns the
    % summarized probability that a contiguous winSize-nt stretch inside the
    % probe binding site is fully unpaired.
    target = upper(strrep(target_seq, '-', '')); target = strrep(target, 'T', 'U');
    probe = upper(strrep(probe_seq, '-', '')); uLen = length(probe);
    winSize = min(winSize, uLen);
    bind_rna_rc = dna_or_rna_to_rna_revcomp(probe);
    hits = strfind(target, bind_rna_rc);
    if isempty(hits), error('Probe reverse complement not found in target.'); end

    start_pos = hits(1); end_pos = start_pos + uLen - 1;
    winStart = max(1, start_pos - targetExt);
    winEnd   = min(length(target), end_pos + targetExt);
    windowU  = target(winStart:winEnd);

    capStart_win = start_pos - winStart + 1;
    capEnd_win   = capStart_win + uLen - 1;
    W = length(windowU); L = max(1, W - 1);

    workdir = tempname; mkdir(workdir);
    fa = fullfile(workdir, 'input.fa');
    fid = fopen(fa, 'w'); fprintf(fid, '>target_window\n%s\n', windowU); fclose(fid);

    rnaPlfold = fullfile(getenv('HOME'), 'viennarna', 'bin', 'RNAplfold');
    cmd = sprintf('cd "%s" && "%s" -W %d -L %d -u %d --temp %.2f --salt %.6f < "%s" > /dev/null 2>&1', ...
        workdir, rnaPlfold, W, L, winSize, tempC, salt_M, fa);
    status = system(cmd);
    if status ~= 0, error('RNAplfold failed.'); end

    d = dir(fullfile(workdir, '*lunp*'));
    if isempty(d), error('No lunp output found.'); end
    lunpFile = fullfile(workdir, d(1).name);

    Pacc = read_lunp_window_stat(lunpFile, capStart_win, capEnd_win, winSize, stat_mode, top_percent);
    R = 8.314; T = tempC + 273.15;
    out = struct('Pacc', Pacc, 'DGopen_kcalmol', (-R * T * log(Pacc)) / 4184);
    rmdir(workdir, 's');
end

function P_stat = read_lunp_window_stat(lunpFile, capStart, capEnd, winSize, stat_mode, top_percent)
    % Parses an RNAplfold _lunp file and summarizes the unpaired probabilities
    % of every winSize-nt stretch fully contained within the binding site.
    lines = readlines(lunpFile);
    hdr_idx = find(contains(lines(1:min(10,end)), "l=1"), 1);

    start_row = hdr_idx + capStart + winSize - 1;
    end_row   = hdr_idx + capEnd;

    P_array = nan(1, max(0, end_row - start_row + 1));
    idx = 1;
    for r = start_row:end_row
        if r > numel(lines), break; end
        tok = regexp(strtrim(lines(r)), '\s+', 'split');
        if numel(tok) > winSize
            v = tok{1 + winSize};
            if ~strcmpi(v, 'NA')
                val = str2double(v);
                if ~isnan(val), P_array(idx) = val; idx = idx + 1; end
            end
        end
    end
    P_array = P_array(1:idx-1);

    if isempty(P_array)
        P_stat = 1e-12;
    else
        if strcmpi(stat_mode, 'avg') || strcmpi(stat_mode, 'mean')
            P_stat = mean(P_array);
        elseif strcmpi(stat_mode, 'max')
            P_stat = max(P_array);
        elseif strcmpi(stat_mode, 'topX') || strcmpi(stat_mode, 'top')
            P_sorted = sort(P_array, 'descend');
            num_top = max(1, ceil(top_percent * length(P_sorted))); % keep at least one
            P_stat = mean(P_sorted(1:num_top));
        elseif strcmpi(stat_mode, 'top25')
            P_sorted = sort(P_array, 'descend');
            P_stat = mean(P_sorted(1:max(1, ceil(0.25 * length(P_sorted)))));
        end
    end
    if isnan(P_stat) || P_stat <= 0, P_stat = 1e-12; end
end

function rc = dna_or_rna_to_rna_revcomp(seq)
    % Reverse complement, returned as RNA.
    s = upper(seq); s = strrep(s, 'U', 'T'); s = strrep(s, 'A', 'x');
    s = strrep(s, 'T', 'A'); s = strrep(s, 'x', 'T'); s = strrep(s, 'C', 'y');
    s = strrep(s, 'G', 'C'); s = strrep(s, 'y', 'G');
    rc = strrep(s(end:-1:1), 'T', 'U');
end
