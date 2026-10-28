%% Accessibility ranking as a function of unpaired-stretch window length
% Generates Supplementary Figure S10.
%
% WHAT IT DOES
%   Scans the HIV-1 consensus genome with a 70-nt sliding window and scores
%   each window by the maximum probability that a contiguous stretch of a
%   given length is fully unpaired anywhere inside it. The scan is repeated
%   for stretch lengths of 10, 6 and 3 nt, and the pol (IN) and env (gp41)
%   capture sites are marked on each panel.
%
%   The point of the figure is that the ranking of the two sites depends on
%   the stretch length chosen. At 10 nt env sits near the bottom of the
%   genome-wide distribution while pol is mid-range; the gap narrows at 6 nt
%   and reverses at 3 nt, where nearly every window in the genome scores
%   above 0.7.
%
%   Note that this is the maximum across the window, not the mean used for
%   the accessibility values reported elsewhere in the manuscript. The
%   maximum is the relevant statistic here because the question is whether a
%   nucleation site of a given length exists anywhere in the binding site,
%   not how open the site is on average.
%
%   RNAplfold is run once in local folding mode over the whole ungapped
%   consensus, which returns the probabilities for every stretch length up to
%   u_max in a single pass. Window scores for all three stretch lengths are
%   then taken from that one result.
%
%   Capture sites are located by searching the consensus for the target
%   sequence of each probe, so no genome coordinates need to be hard-coded.
%
% REQUIREMENTS
%   MATLAB R2020b or later (uses readlines, movmax and exportgraphics).
%   ViennaRNA 2.6.4 or later, installed locally. RNAplfold is called as an
%   external executable, so the ViennaRNA binaries must be present on disk.
%   Download: https://www.tbi.univie.ac.at/RNA/
%
%   Input file: Los_HIV_Full_consensus.txt, the same FASTA consensus used by
%   P_profile_T_Sweep.m.
%
% SETUP
%   By default the script expects RNAplfold at:
%       ~/viennarna/bin/RNAplfold
%   If ViennaRNA is installed elsewhere, edit the rnaPlfold path inside the
%   function run_rnaplfold_local near the bottom of this file.
%   On Windows, point it to RNAplfold.exe.
%
% HOW TO RUN
%   Put Los_HIV_Full_consensus.txt in the working directory, open this file
%   in MATLAB and press Run, or call it by name from the command window.
%   No input arguments are required; all settings live in section 1 below.
%   RNAplfold is called once, so the whole script finishes in a few minutes.
%   Outputs: one PNG at 300 dpi, plus a per-site summary printed to the
%   command window.
%
% REPRODUCIBILITY
%   The random seed is fixed below so that the IUPAC ambiguity substitutions
%   are identical on repeated runs.

clc;
clear;
close all;

rng(0);  % fixed seed for the IUPAC ambiguity substitutions

outputPNG = 'Window_Length_Accessibility_Scan.png';

%% ===================== 1. INPUT & SETTINGS =====================
consensusFile = 'Los_HIV_Full_consensus.txt';

windowLen     = 70;             % capture window length (nt)
segLens       = [10, 6, 3];     % unpaired stretch lengths to compare
u_max         = max(segLens);   % RNAplfold computes every length up to this
tempC         = 65;             % folding temperature (C)
sodiumC       = 0.5;            % [Na+] (M)

% Local folding span. RNAplfold restricts base pairs to within localSpan
% nucleotides, so each position effectively sees this much sequence context
% on either side. Set to match the flanking length quoted in the figure
% caption.
localSpan     = 200;

maxMismatches = 5;              % tolerance when locating capture sites

% Capture sites, given as the target sequence on the viral RNA, that is the
% reverse complement of the capture probe listed in Supplementary Table S1
site_seqs = {
    'GCAGTACAAATGGCAGTATTCATTCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGG', ... % pol (IN)
    'TGGGAGCAGCAGGAAGCACTATGGGCGCGGCGTCAATAACGCTGACGGTACAGGCCAGACAATTATTGTC'  ... % env (gp41)
};
site_labels = {'pol (IN)', 'env (gp41)'};
site_colors = [0.93 0.49 0.19;   % orange
               0.85 0.11 0.11];  % red

%% ===================== 2. LOAD AND PREPARE SEQUENCE =====================
sequences = callSeq(consensusFile);
consensusRaw = sequences{1};

% Remove alignment gaps, then resolve IUPAC ambiguity codes
consensusNoGap = strrep(consensusRaw, '-', '');
consensus = resolveDegenerateOne(consensusNoGap);
L_total = length(consensus);

fprintf('Consensus loaded: %d nt after gap removal\n', L_total);
fprintf('Folding at %.0f C, [Na+] = %.2f M, local span %d nt\n\n', tempC, sodiumC, localSpan);

%% ===================== 3. FOLD ONCE, SCORE EVERY WINDOW LENGTH =====================
% P(i, u) is the probability that the u-nt stretch ENDING at position i is
% fully unpaired. One RNAplfold call returns every u from 1 to u_max.
fprintf('Running RNAplfold over the full consensus...\n');
P = run_rnaplfold_local(consensus, tempC, sodiumC, localSpan, u_max);
fprintf('Done.\n\n');

numWindows = L_total - windowLen + 1;
window_score = nan(numWindows, numel(segLens));

for s = 1:numel(segLens)
    u = segLens(s);

    % Convert to "stretch STARTING at i": the stretch starting at i ends at
    % i + u - 1, so p_start(i) = P(i + u - 1, u)
    p_start = nan(L_total, 1);
    p_start(1:L_total-u+1) = P(u:L_total, u);

    % Maximum over every stretch start that falls fully inside the 70-nt
    % window beginning at j, that is i from j to j + windowLen - u
    n_sub = windowLen - u + 1;
    rolling = movmax(p_start, [0 n_sub-1], 'omitnan');
    window_score(:, s) = rolling(1:numWindows);
end

%% ===================== 4. LOCATE THE CAPTURE SITES =====================
site_start = nan(numel(site_seqs), 1);

fprintf('=== Capture site positions (maximum mismatches allowed: %d) ===\n', maxMismatches);
for k = 1:numel(site_seqs)
    query = upper(strrep(site_seqs{k}, '-', ''));
    [pos, mm] = find_best_match(consensus, query);

    if mm <= maxMismatches
        site_start(k) = pos;
        fprintf('[%s] %d-%d nt, %d mismatches\n', ...
            site_labels{k}, pos, pos + length(query) - 1, mm);
    else
        warning('[%s] no match found (best Hamming distance was %d)', site_labels{k}, mm);
    end
end
fprintf('\n');

%% ===================== 5. FIGURE =====================
fig = figure('Color', 'w', 'Name', 'Accessibility vs window length', ...
             'Position', [100, 100, 1000, 820]);
tl = tiledlayout(numel(segLens), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:numel(segLens)
    u = segLens(s);
    ax = nexttile;
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

    plot(1:numWindows, window_score(:, s), '-', 'Color', [0.30 0.55 0.75], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Mark each capture site at the window that starts exactly on it.
    % Markers only, with no labels or annotation, so that the panels can be
    % labeled by hand afterwards.
    for k = 1:numel(site_seqs)
        if isnan(site_start(k)), continue; end
        j = site_start(k);
        if j < 1 || j > numWindows, continue; end
        plot(j, window_score(j, s), 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', site_colors(k,:), 'MarkerEdgeColor', 'none');
    end

    ylim([0 1]);
    xlim([1 numWindows]);
    ylabel(sprintf('Max P(unpaired), %d nt', u), 'FontSize', 14);
    set(ax, 'FontSize', 14);

    if s < numel(segLens)
        set(ax, 'XTickLabel', []);
    end
end

xlabel(tl, 'Window start position (consensus, 1-based, ungapped)', 'FontSize', 14);

exportgraphics(fig, outputPNG, 'Resolution', 300);
fprintf('Figure saved to %s (300 dpi)\n\n', outputPNG);

%% ===================== 6. SITE SUMMARY =====================
fprintf('=== Capture site scores and genome-wide percentile rank ===\n');
fprintf('%-12s', 'Site');
for s = 1:numel(segLens)
    fprintf('%18s', sprintf('%d-nt', segLens(s)));
end
fprintf('\n');

for k = 1:numel(site_seqs)
    if isnan(site_start(k)), continue; end
    j = site_start(k);
    fprintf('%-12s', site_labels{k});
    for s = 1:numel(segLens)
        v = window_score(j, s);
        pct = mean(window_score(:, s) <= v, 'omitnan') * 100;
        fprintf('%10.3f (%4.1f)', v, pct);
    end
    fprintf('\n');
end
fprintf('Values are max P(unpaired); percentile rank genome-wide in parentheses.\n');

for s = 1:numel(segLens)
    frac = mean(window_score(:, s) > 0.7, 'omitnan') * 100;
    fprintf('At %2d nt, %.1f%% of genome-wide windows exceed 0.7\n', segLens(s), frac);
end
fprintf('\n');

%% ===================== 7. HELPER FUNCTIONS =====================
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
    iupacMap = containers.Map( ...
        {'A','C','G','T','U','R','Y','S','W','K','M','B','D','H','V','N','?'}, ...
        {'A','C','G','T','T','AG','CT','CG','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT','ACGT'});
    seq = repmat('A', size(degSeq));
    for i = 1:length(degSeq)
        code = upper(degSeq(i));
        if ~iupacMap.isKey(code), error('Invalid IUPAC code "%s"', code); end
        letters = iupacMap(code);
        seq(i) = letters(randi(numel(letters)));
    end
end

function [bestPos, bestMM] = find_best_match(target, query)
    % Sliding-window search for the position with the smallest Hamming
    % distance to the query.
    t = upper(strrep(target, 'U', 'T'));
    q = upper(strrep(query, 'U', 'T'));
    L = length(q);
    bestMM = Inf; bestPos = -1;
    for i = 1:(length(t) - L + 1)
        mm = sum(t(i:i+L-1) ~= q);
        if mm < bestMM
            bestMM = mm; bestPos = i;
            if mm == 0, break; end
        end
    end
end

function P = run_rnaplfold_local(seq, tempC, salt_M, localSpan, u_max)
    % Runs RNAplfold once in local folding mode over the whole sequence and
    % returns the unpaired-probability matrix.
    % P(i, u) is the probability that the u-nt stretch ending at position i
    % is fully unpaired. Entries that do not exist are left as NaN.
    s = upper(strrep(seq, 'T', 'U'));

    workdir = tempname; mkdir(workdir);
    cleanupObj = onCleanup(@() rmdir(workdir, 's'));

    fa = fullfile(workdir, 'input.fa');
    fid = fopen(fa, 'w'); fprintf(fid, '>consensus\n%s\n', s); fclose(fid);

    rnaPlfold = fullfile(getenv('HOME'), 'viennarna', 'bin', 'RNAplfold');
    cmd = sprintf('cd "%s" && "%s" -W %d -L %d -u %d --temp %.2f --salt %.6f < "%s" > /dev/null 2>&1', ...
        workdir, rnaPlfold, localSpan, localSpan, u_max, tempC, salt_M, fa);
    status = system(cmd);
    if status ~= 0, error('RNAplfold failed.'); end

    d = dir(fullfile(workdir, '*lunp*'));
    if isempty(d), error('No lunp output found.'); end
    P = read_lunp_matrix(fullfile(workdir, d(1).name), length(s), u_max);
end

function P = read_lunp_matrix(lunpFile, L_total, u_max)
    % Parses a complete RNAplfold _lunp file into an L_total by u_max matrix.
    lines = readlines(lunpFile);
    hdr_idx = find(contains(lines(1:min(10, end)), "l=1"), 1);
    if isempty(hdr_idx), error('Could not find the lunp header row.'); end

    P = nan(L_total, u_max);
    for r = (hdr_idx + 1):numel(lines)
        tok = regexp(strtrim(lines(r)), '\s+', 'split');
        if numel(tok) < 2, continue; end
        i = str2double(tok{1});
        if isnan(i) || i < 1 || i > L_total, continue; end
        for u = 1:min(u_max, numel(tok) - 1)
            v = tok{1 + u};
            if ~strcmpi(v, 'NA')
                val = str2double(v);
                if ~isnan(val), P(i, u) = val; end
            end
        end
    end
end
