%% Genome-wide accessibility landscape of the HIV-1 consensus sequence
% Step 1 of 2. Computes the data underlying Supplementary Figure S4C.
% Run Saved_Data_to_Figure.m afterwards to render the published figure.
%
% WHAT IT DOES
%   Scans the full-length HIV-1 consensus sequence with a sliding window,
%   shifted in 10-nt increments. At each position the window is treated as a
%   candidate capture probe binding site and is extracted together with
%   150-nt flanking sequence on both sides, so that local structural context
%   is preserved. RNAplfold is then used to compute the probability that a
%   contiguous 10-nt stretch inside the binding site is fully unpaired, and
%   these probabilities are averaged across the site to give its accessibility
%   (P_unpaired). The scan is repeated at every temperature from 25 to 85 C in
%   5 C steps, and for four candidate probe lengths (19, 40, 70 and 120 nt).
%
%   Alignment gaps are removed before folding. IUPAC ambiguity codes and
%   unresolved bases are each replaced at random by one of the canonical bases
%   they represent. Windows in which at least half the bases are undefined are
%   skipped and left as NaN.
%
%   Results are written to HIV_Thermal_Data_Saved2.mat. The plotting script
%   reads that file, so the long computation only has to be run once.
%
% REQUIREMENTS
%   MATLAB R2020b or later (uses readlines and exportgraphics).
%   ViennaRNA 2.6.4 or later, installed locally. RNAplfold is called as an
%   external executable, so the ViennaRNA binaries must be present on disk.
%   Download: https://www.tbi.univie.ac.at/RNA/
%
%   Input file: Los_HIV_Full_consensus.txt, a FASTA file holding the HIV-1
%   consensus sequence. The published analysis used the Los Alamos National
%   Laboratory HIV sequence database (HIV-1/SIVcpz, Region: GENOME,
%   Subtype: M group with CRFs, 2021), from which the consensus was generated
%   by taking the majority nucleotide at each alignment position.
%
% SETUP
%   By default the script expects RNAplfold at:
%       ~/viennarna/bin/RNAplfold
%   If ViennaRNA is installed elsewhere, edit the rnaPlfold path inside the
%   function DGopen_stat_for_probe near the bottom of this file.
%   On Windows, point it to RNAplfold.exe.
%
% HOW TO RUN
%   Put Los_HIV_Full_consensus.txt in the working directory, open this file
%   in MATLAB and press Run, or call it by name from the command window.
%   No input arguments are required; all settings live in section 1 below.
%   The full scan makes roughly four million RNAplfold calls and takes about
%   a day on a desktop machine. To test the pipeline quickly, set
%   test_limit_N to a few thousand so that only the start of the genome is
%   scanned.
%   Outputs: HIV_Thermal_Data_Saved2.mat plus four preview PNGs at 300 dpi,
%   named in the Filename1-4 variables below. The published figure is
%   produced by Saved_Data_to_Figure.m, not by the previews here.
%
% REPRODUCIBILITY
%   The random seed is fixed below so that the ambiguity-code substitutions
%   are identical on repeated runs. Note that the published .mat file was
%   generated before the seed was fixed, so a fresh run will differ from it
%   very slightly at ambiguous positions.

clc;
clear;
close all;

rng(0);  % fixed seed for the IUPAC ambiguity substitutions

Filename1 = 'P_profile_Across_avg_10nt_HIV1.png';
Filename2 = 'P_profile_Across_avg_10nt_HIV2.png';
Filename3 = 'P_profile_Across_avg_10nt_HIV3.png';
Filename4 = 'P_profile_Across_avg_10nt_HIV4.png';

%% ===================== 1. INPUT & SETTINGS =====================
consensusFile = 'Los_HIV_Full_consensus.txt';

% Main parameters
test_limit_N = inf;  % number of positions to scan; set to a small value to test
step_size = 10;      % sliding window shift (nt)
probeLength = 70;    % overwritten in the loop by probeLength_range below
unpair_window = 10;  % length of the contiguous unpaired stretch scored

% Statistic used to summarize P_unpaired across the binding site
stat_mode = 'avg';   % 'avg', 'max', or 'topX'
top_percent = 0.10;  % used only when stat_mode is 'topX'
sodiumC = 0.5;       % [Na+] (M), matched to the experimental conditions
extension = 150;     % flanking sequence retained on each side (nt)

% Temperature sweep
T_min = 25;
T_max = 85;
dT_calc = 5;
T_sweep = T_min:dT_calc:T_max;
num_T_calc = length(T_sweep);

% Temperatures drawn in the landscape preview (Figure 1)
T_show = [25, 65];
num_T_show = length(T_show);

%% ===================== 2. LOAD SEQUENCE =====================
sequences = callSeq(consensusFile);
consensusSequence = sequences{1};
seqLength = length(consensusSequence);
actual_limit = min(test_limit_N, seqLength);

% Candidate probe lengths analysed in parallel
probeLength_range = [19, 40, 70, 120];
num_cases = length(probeLength_range);

% Window start positions to scan
scan_indices = 1:step_size:actual_limit;
num_scan_points = length(scan_indices);

% Total step count, used only for the progress readout
total_steps = num_cases * num_scan_points;
current_step = 0;

fprintf('Genome-wide accessibility landscape and global temperature trend\n');
fprintf('Probe lengths analysed: '); disp(probeLength_range);
fprintf('Step: %d nt, temperature range: %d-%d C\n', step_size, T_min, T_max);
fprintf('Statistic: %s (mean of top %d%%)\n', stat_mode, round(top_percent * 100));
fprintf('Total steps: %d\n\n', total_steps);

% Container for the results, one matrix per probe length
P_landscape_all = cell(num_cases, 1);

%% ===================== 3. DATA COLLECTION LOOP =====================
for idxLen = 1:num_cases
    probeLength = probeLength_range(idxLen);
    fprintf('=== [%d nt probe] starting data collection ===\n', probeLength);

    P_landscape_temp = nan(num_T_calc, num_scan_points);

    for pt_idx = 1:num_scan_points
        i = scan_indices(pt_idx);
        current_step = current_step + 1;

        % Walk forward until probeLength non-gap bases have been collected
        valid_bases = 0; endIndex = i - 1;
        while valid_bases < probeLength && endIndex < seqLength
            endIndex = endIndex + 1;
            if consensusSequence(endIndex) ~= '-', valid_bases = valid_bases + 1; end
        end
        if valid_bases < probeLength
            % Ran off the end of the sequence; stop scanning this probe length
            current_step = current_step + (num_scan_points - pt_idx);
            break;
        end

        cap_region_raw = consensusSequence(i:endIndex);
        cap_region_nogap = strrep(cap_region_raw, '-', '');

        % Skip windows that are at least half undefined
        cap_upper = upper(cap_region_nogap);
        num_unknown = sum(cap_upper == 'N' | cap_upper == '?');
        if num_unknown / probeLength >= 0.5, continue; end

        % Extract the binding site together with its flanking context
        start_ext = max(1, i - extension);
        end_ext   = min(seqLength, endIndex + extension);
        extended_raw = consensusSequence(start_ext:end_ext);
        extended_nogap = strrep(extended_raw, '-', '');

        % Resolve IUPAC ambiguity codes to canonical bases
        extended_resolved = resolveDegenerateOne(extended_nogap);

        % Locate the binding site within the gap-free extended sequence
        left_flank_raw = consensusSequence(start_ext : i - 1);
        left_flank_nogap = strrep(left_flank_raw, '-', '');
        cap_start_idx = length(left_flank_nogap) + 1;

        cap_region_resolved = extended_resolved(cap_start_idx : cap_start_idx + probeLength - 1);
        probe_rna = strrep(strrep(upper(reverseComplementIUPAC(cap_region_resolved)), 'U', 'T'), 'T', 'U');

        try
            for t_idx = 1:num_T_calc
                tempC = T_sweep(t_idx);
                out_v = DGopen_stat_for_probe(extended_resolved, probe_rna, tempC, sodiumC, extension, unpair_window, stat_mode, top_percent);
                P_landscape_temp(t_idx, pt_idx) = out_v.Pacc;
            end
        catch ME
            disp(['Error at position ', num2str(i), ': ', ME.message]);
        end

        if mod(pt_idx, 10) == 0 || pt_idx == num_scan_points
            overall_percent = (current_step / total_steps) * 100;
            fprintf('  > Pos %d / %d scanned ... (overall progress: %5.1f%%)\n', i, actual_limit, overall_percent);
        end
    end
    P_landscape_all{idxLen} = P_landscape_temp;
    fprintf(' [%d nt probe] complete.\n\n', probeLength);
end

fprintf('Data collection complete. Saving results...\n');

% =========================================================================
% Save everything the plotting script needs, so the scan is run only once
% =========================================================================
save('HIV_Thermal_Data_Saved2.mat', 'P_landscape_all', 'T_sweep', ...
    'scan_indices', 'actual_limit', 'T_show', 'num_T_show', ...
    'probeLength_range', 'num_cases', 'T_min', 'T_max');
fprintf('Saved: HIV_Thermal_Data_Saved2.mat\n\n');


%% ===================== 4. PREVIEW PLOTS =====================
% Quick previews of the data just computed. The published figure is rendered
% by Saved_Data_to_Figure.m from the saved .mat file.
trend_colors = lines(4);
idx_70 = find(probeLength_range == 70, 1);
idx_19 = find(probeLength_range == 19, 1);

% -----------------------------------------------------------------
% [Figure 1] Genome-wide landscape (70 nt, two representative temperatures)
% -----------------------------------------------------------------
colors_show = jet(num_T_show);
fig1 = figure('Name', 'Fig 1: Landscape (70nt)');
box on; grid on; hold on;
P_land_70 = P_landscape_all{idx_70};

for t_idx = 1:num_T_show
    real_idx = find(T_sweep == T_show(t_idx));
    if ~isempty(real_idx)
        plot(scan_indices, P_land_70(real_idx, :), '-o', 'Color', colors_show(t_idx, :), ...
            'LineWidth', 2.5, 'MarkerSize', 4, 'DisplayName', sprintf('%d \\circC', T_show(t_idx)));
    end
end
set(gca, 'FontSize', 16);
xlabel('Position on HIV Sequence (nt)');
ylabel('Unpaired Prob');
legend('Location', 'northwest', 'FontSize', 14);
ylim([-0.05 1.05]);
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
legend('Location', 'northwest', 'FontSize', 14);
xticks(T_min:10:T_max); xlim([25 85]);
exportgraphics(gcf, Filename2, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 3] Global temperature trend, 19 nt vs 70 nt
% -----------------------------------------------------------------
fig3 = figure('Name', 'Fig 3: Trend (19nt vs 70nt)');
box on; grid on; hold on;
compare_idx = [idx_19, idx_70];
comp_colors = [trend_colors(1,:); trend_colors(3,:)];

% Shaded SD bands first, mean lines on top
for k = 1:2
    c_idx = compare_idx(k);
    P_mat = P_landscape_all{c_idx};
    g_mean = mean(P_mat, 2, 'omitnan')';
    g_std  = std(P_mat, 0, 2, 'omitnan')';
    u_b = min(1, g_mean + g_std);
    l_b = max(0, g_mean - g_std);
    fill([T_sweep, fliplr(T_sweep)], [u_b, fliplr(l_b)], comp_colors(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
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
xticks(T_min:10:T_max); xlim([25 85]);
exportgraphics(gcf, Filename3, 'Resolution', 300);

% -----------------------------------------------------------------
% [Figure 4] Global temperature trend, all four probe lengths
% -----------------------------------------------------------------
fig4 = figure('Name', 'Fig 4: Trend (All Lengths)');
box on; grid on; hold on;

% Shaded SD bands first, mean lines on top
for k = 1:num_cases
    P_mat = P_landscape_all{k};
    g_mean = mean(P_mat, 2, 'omitnan')';
    g_std  = std(P_mat, 0, 2, 'omitnan')';
    u_b = min(1, g_mean + g_std);
    l_b = max(0, g_mean - g_std);
    fill([T_sweep, fliplr(T_sweep)], [u_b, fliplr(l_b)], trend_colors(k,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
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
xticks(T_min:10:T_max); xlim([25 85]);
exportgraphics(gcf, Filename4, 'Resolution', 300);

fprintf('All preview images (Fig 1-4) exported.\n');

%% ===================== 5. HELPER FUNCTIONS =====================
function sequences = callSeq(filename)
    % Minimal FASTA reader. Returns a cell array of sequence strings.
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open the file: %s', filename); end
    sequences = {}; currentSequence = '';
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line), continue; end
        if line(1) == '>'
            if ~isempty(currentSequence), sequences{end + 1} = currentSequence; currentSequence = ''; end
        else
            currentSequence = [currentSequence line];
        end
    end
    if ~isempty(currentSequence), sequences{end + 1} = currentSequence; end
    fclose(fid);
end

function seq = resolveDegenerateOne(degSeq)
    % Replaces each IUPAC ambiguity code with one of the canonical bases it
    % represents, chosen at random and independently at every position.
    % Gaps are left in place; they are stripped by the caller.
    iupacMap = containers.Map( ...
        {'A','C','G','T','R','Y','S','W','K','M','B','D','H','V','N','-','?'}, ...
        {'A','C','G','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT','-','ACGT'});
    seq = repmat('-', size(degSeq));
    for i = 1:length(degSeq)
        code = upper(degSeq(i));
        if ~iupacMap.isKey(code), error('Invalid IUPAC code "%s"', code); end
        letters = iupacMap(code);
        seq(i) = letters(randi(numel(letters)));
    end
end

function rc = reverseComplementIUPAC(seq)
    % Reverse complement, returned as DNA.
    s = upper(seq); s = strrep(s, 'U', 'T'); s = strrep(s, 'A', 'x');
    s = strrep(s, 'T', 'A'); s = strrep(s, 'x', 'T'); s = strrep(s, 'C', 'y');
    s = strrep(s, 'G', 'C'); s = strrep(s, 'y', 'G'); rc = s(end:-1:1);
end

function out = DGopen_stat_for_probe(target_seq, probe_seq, tempC, salt_M, targetExt, winSize, stat_mode, top_percent)
    % Runs RNAplfold on the supplied target sequence and returns the
    % summarized probability that a contiguous winSize-nt stretch inside the
    % probe binding site is fully unpaired.
    % Note: targetExt is accepted for interface compatibility but is not used
    % here, because the caller has already added the flanking context.
    target = upper(strrep(target_seq, '-', '')); target = strrep(target, 'T', 'U');
    probe = upper(strrep(probe_seq, '-', '')); uLen = length(probe);
    winSize = min(winSize, uLen);
    bind_rna_rc = strrep(upper(reverseComplementIUPAC(probe)), 'T', 'U');
    hits = strfind(target, bind_rna_rc);
    if isempty(hits), error('Probe not found.'); end

    start_pos = hits(1); end_pos = start_pos + uLen - 1;
    W = length(target); L = max(1, W - 1);

    workdir = tempname; mkdir(workdir);
    fa = fullfile(workdir, 'input.fa');
    fid = fopen(fa, 'w'); fprintf(fid, '>target_window\n%s\n', target); fclose(fid);

    rnaPlfold = fullfile(getenv('HOME'), 'viennarna', 'bin', 'RNAplfold');
    cmd = sprintf('cd "%s" && "%s" -W %d -L %d -u %d --temp %.2f --salt %.6f < "%s" > /dev/null 2>&1', ...
        workdir, rnaPlfold, W, L, winSize, tempC, salt_M, fa);
    status = system(cmd);
    if status ~= 0, error('RNAplfold failed.'); end
    d = dir(fullfile(workdir, '*lunp*')); lunpFile = fullfile(workdir, d(1).name);

    P_array = read_lunp_window_array(lunpFile, start_pos, end_pos, winSize);
    if strcmpi(stat_mode, 'avg') || strcmpi(stat_mode, 'mean')
        Pacc = mean(P_array);
    elseif strcmpi(stat_mode, 'max')
        Pacc = max(P_array);
    elseif strcmpi(stat_mode, 'topX') || strcmpi(stat_mode, 'top')
        P_sorted = sort(P_array, 'descend');
        num_top = max(1, ceil(top_percent * length(P_sorted)));
        Pacc = mean(P_sorted(1:num_top));
    else
        Pacc = mean(P_array);
    end

    out = struct('Pacc', Pacc);
    rmdir(workdir, 's');
end

function P_array = read_lunp_window_array(lunpFile, capStart, capEnd, winSize)
    % Parses an RNAplfold _lunp file and returns the unpaired probabilities of
    % every winSize-nt stretch fully contained within the binding site.
    lines = readlines(lunpFile);
    hdr_idx = find(contains(lines(1:min(10,end)), "l=1"), 1);
    start_row = hdr_idx + capStart + winSize - 1; end_row = hdr_idx + capEnd;

    P_array = nan(1, max(0, end_row - start_row + 1)); idx = 1;
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
    if isempty(P_array), P_array = 1e-12; end
end
