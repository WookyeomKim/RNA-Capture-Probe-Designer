%% Experimental structure probing of the HIV-1 genome (SHAPE reactivity)
% Generates Supplementary Figure S9.
%
% WHAT IT DOES
%   Loads a published per-nucleotide SHAPE reactivity profile for the HIV-1
%   genome and locates each capture probe binding site within it by
%   sliding-window sequence matching, keeping the position with the smallest
%   Hamming distance. For context the profile is drawn as a 75-nt moving
%   median, chosen to match the scale of the 70-nt capture footprint.
%
%   Each capture site is then marked on that smoothed profile. Every position
%   inside the 70-nt window is drawn individually, with a shaded band spanning
%   the values the profile takes across the window, so that the stretch of the
%   genome-wide curve occupied by the site is visible rather than collapsed
%   into a single point. Raw reactivity is deliberately not plotted, because
%   individual nucleotides scatter over the full 0 to 2 range and would hide
%   where the site sits; the unsmoothed mean, median, SD and range are
%   reported in the console instead.
%
%   No summary line is drawn across the window. The smoothed profile and the
%   unsmoothed median are different quantities, and drawing one as a level to
%   be read against the percentile grid would invite comparison with the
%   other. The figure therefore shows position only, and the manuscript
%   quotes the percentile rank computed below.
%
%   The unsmoothed median is additionally converted into a genome-wide
%   percentile rank, computed by taking the median across every 70-nt window
%   in the genome. The percentile rank is the value quoted in the manuscript.
%
%   Reactivities above 2.0 or below 0 are treated as outliers and excluded.
%
% REQUIREMENTS
%   MATLAB R2020b or later (uses exportgraphics and yline/xline labels).
%   The Statistics and Machine Learning Toolbox is needed for prctile.
%
%   Input file: HIV1_NL43_SHAPE.txt, a tab-delimited table with four columns
%   and no header: position, nucleotide, SHAPE reactivity, pairing
%   probability. The published analysis used the HIV-1 NL4-3 genome SHAPE
%   dataset of Watts et al., Nature 460, 711-716 (2009).
%
% HOW TO RUN
%   Put HIV1_NL43_SHAPE.txt in the working directory, open this file in
%   MATLAB and press Run, or call it by name from the command window.
%   No input arguments are required. Capture site sequences are listed in the
%   seqs_raw cell array in section 2; edit that array to evaluate different
%   sites. Gaps written as "-" are stripped before matching.
%   Runtime is well under a minute.
%   Outputs: one PNG at 300 dpi, plus a per-target summary printed to the
%   command window.
%
% NOTE ON PERCENTILE DIRECTION
%   Two complementary numbers are reported for each site. The percentile rank
%   is the fraction of genome-wide windows that are no more reactive than this
%   one, so a low value means an unusually unreactive, highly structured site.
%   The "top" figure is its complement, used for the figure legend. The
%   manuscript quotes the percentile rank.

clc;
clear;
close all;

%% ===================== 1. LOAD AND PARSE DATA =====================
filename = 'HIV1_NL43_SHAPE.txt';

opts = delimitedTextImportOptions("NumVariables", 4);
opts.Delimiter = "\t";
opts.VariableNames = ["Position", "Nucleotide", "SHAPE_Reactivity", "Pairing_Probability"];
opts.VariableTypes = ["double", "string", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

data = readtable(filename, opts);

% Treat reactivities above 2.0 and negative values as outliers
data.SHAPE_Reactivity(data.SHAPE_Reactivity > 2.0 | data.SHAPE_Reactivity < 0) = NaN;

fullGenomeSeq = strjoin(data.Nucleotide, '');
genomeLen = strlength(fullGenomeSeq);

%% ===================== 2. TARGET SEQUENCES AND SEARCH SETTINGS =====================
maxMismatches = 5;   % maximum Hamming distance accepted as a match

seqs_raw = {
    'TGGGAGCAGCAGGAAGCACTATGGGCGCGGCGTCAATAACGCTGACGGTACAGGCCAGACAATTATTGTC', ... % Target 1 (Red)
    'AGCAGTACAAATGGCAGTATTCATTCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGG', ... % Target 2 (Orange)
    'AAATATAGACTAAA-ACATTTAGTA----TGGGCAAGCAGGGAGCTGGAAAGAT-TCGCACTTAACCCTGGCCTTT' ... % Target 3 (Green)
};

labels = {'Target 1 (Red)', 'Target 2 (Orange)', 'Target 3 (Green)'};
colors = [1 0 0; 1 0.5 0; 0 0.7 0]; % red, orange, green

% Preallocate so that a failed match leaves a clearly empty entry
matches = struct('start', {}, 'end', {}, 'mid', {}, 'mismatches', {}, 'label', {}, 'color', {});

fprintf('=== Sequence matching (maximum mismatches allowed: %d) ===\n', maxMismatches);

for k = 1:length(seqs_raw)
    cleanSeq = upper(strrep(strrep(strtrim(seqs_raw{k}), '-', ''), 'T', 'U'));
    L = strlength(cleanSeq);

    minMismatchesFound = Inf;
    bestStartPos = -1;

    % Slide the query across the genome, keeping the best-matching position
    for i = 1:(genomeLen - L + 1)
        windowSeq = extractBetween(fullGenomeSeq, i, i + L - 1);
        mismatches = sum(char(windowSeq) ~= char(cleanSeq));

        if mismatches < minMismatchesFound
            minMismatchesFound = mismatches;
            bestStartPos = i;
        end
    end

    if minMismatchesFound <= maxMismatches
        idx = numel(matches) + 1;
        matches(idx).start = bestStartPos;
        matches(idx).end = bestStartPos + L - 1;
        matches(idx).mid = round((bestStartPos + matches(idx).end) / 2);
        matches(idx).mismatches = minMismatchesFound;
        matches(idx).label = labels{k};
        matches(idx).color = colors(k, :);

        fprintf('[%s] matched at %d-%d nt (centre %d nt), %d mismatches\n', ...
            labels{k}, matches(idx).start, matches(idx).end, matches(idx).mid, minMismatchesFound);
    else
        warning('[%s] no match found (best Hamming distance was %d)', labels{k}, minMismatchesFound);
    end
end
fprintf('\n');

% 75-nt moving median, matched to the scale of the 70-nt capture footprint
smoothSHAPE = movmedian(data.SHAPE_Reactivity, 75, 'omitnan');

%% ===================== 3. GENOME-WIDE 70-nt WINDOW DISTRIBUTION =====================
% Percentiles are computed from unsmoothed reactivity, taking the median
% across every 70-nt window in the genome.
winLen = 70;
numWindows = genomeLen - winLen + 1;
genome_70nt_medians = zeros(numWindows, 1);

for i = 1:numWindows
    genome_70nt_medians(i) = median(data.SHAPE_Reactivity(i : i + winLen - 1), 'omitnan');
end

% Reference levels for the horizontal grid lines in the figure.
% Higher SHAPE reactivity corresponds to a more accessible window.
pct_levels = [99, 90, 80, 70, 60, 50, 40, 30, 20, 10];
pct_labels = {'Top 1%', 'Top 10%', 'Top 20%', 'Top 30%', 'Top 40%', ...
              'Top 50%', 'Top 60%', 'Top 70%', 'Top 80%', 'Top 90%'};

pct_values = prctile(genome_70nt_medians, pct_levels);

fprintf('=== Genome-wide 70-nt median SHAPE percentile reference levels ===\n');
for p = 1:length(pct_levels)
    fprintf('  - %-8s accessibility threshold: SHAPE >= %.3f\n', pct_labels{p}, pct_values(p));
end
fprintf('\n');

%% ===================== 4. FIGURE: PERCENTILE GRID AND TARGET WINDOWS =====================
fig1 = figure('Name', 'HIV-1 Genome SHAPE Profile with Percentile Thresholds', ...
              'NumberTitle', 'off', 'Position', [100, 100, 1300, 650]);

% 1) 75-nt smoothed profile
hPlot = plot(data.Position, smoothSHAPE, 'Color', [0.2 0.5 0.8], 'LineWidth', 1.2);
hold on;

% 2) Horizontal percentile reference lines
grayColor = [0.6 0.6 0.6];
for p = 1:length(pct_levels)
    if strcmp(pct_labels{p}, 'Top 1%')
        % Top 1% emphasised with a purple dotted line
        yline(pct_values(p), ':', 'Top 1%', 'Color', [0.5 0 0.5], 'LineWidth', 1.2, ...
            'LabelHorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    elseif strcmp(pct_labels{p}, 'Top 50%')
        % Top 50% (the genome-wide median) emphasised with a black dashed line
        yline(pct_values(p), '--', 'Top 50% (Median)', 'Color', [0 0 0], 'LineWidth', 1.0, ...
            'LabelHorizontalAlignment', 'right', 'FontSize', 9);
    else
        % Remaining levels drawn as faint gray dotted lines
        yline(pct_values(p), ':', pct_labels{p}, 'Color', grayColor, 'LineWidth', 0.8, ...
            'LabelHorizontalAlignment', 'right', 'FontSize', 8);
    end
end

% 3) Capture windows: every position shown individually, with the range and
%    mean of the window drawn on top
legendHandles = [hPlot];
legendLabels = {'75-nt Smoothed SHAPE'};

for k = 1:numel(matches)
    m = matches(k);

    % Window boundaries
    xline(m.start, ':', 'Color', [m.color 0.6], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(m.end, ':', 'Color', [m.color 0.6], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % Smoothed profile across the window. The capture site is marked on the
    % 75-nt moving median rather than on raw reactivity, because individual
    % nucleotides scatter over the full 0 to 2 range and would obscure where
    % the site sits relative to the genome-wide profile.
    target_pos    = data.Position(m.start:m.end);
    target_smooth = smoothSHAPE(m.start:m.end);
    valid = ~isnan(target_smooth);

    smooth_min  = min(target_smooth, [], 'omitnan');
    smooth_max  = max(target_smooth, [], 'omitnan');

    % Unsmoothed statistics for the window, reported but not plotted
    target_raw    = data.SHAPE_Reactivity(m.start:m.end);
    target_mean   = mean(target_raw, 'omitnan');
    target_median = median(target_raw, 'omitnan');
    target_min    = min(target_raw, [], 'omitnan');
    target_max    = max(target_raw, [], 'omitnan');
    target_sd     = std(target_raw, 'omitnan');
    n_valid       = sum(~isnan(target_raw));

    % Shaded band spanning the smoothed values across the window
    patch([m.start m.end m.end m.start], ...
          [smooth_min smooth_min smooth_max smooth_max], ...
          m.color, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % One marker per nucleotide position, lying on the smoothed curve.
    % No summary line is drawn across the window: the points show where the
    % site sits on the profile, while the numerical summary quoted in the
    % manuscript is the unsmoothed median, reported in the console below.
    hSite = scatter(target_pos(valid), target_smooth(valid), 16, m.color, 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');

    % Genome-wide percentile rank of this window, based on its median.
    % Low rank means an unusually unreactive, highly structured site.
    pct_rank = mean(genome_70nt_medians <= target_median, 'omitnan') * 100;
    top_pct  = 100 - pct_rank;

    legendHandles(end+1) = hSite;
    legendLabels{end+1} = sprintf('%s (%d-%d nt, %.1fth percentile)', ...
        m.label, m.start, m.end, pct_rank);

    % Console summary
    fprintf('[%s] position: %d-%d nt (%d of %d positions with data)\n', ...
        m.label, m.start, m.end, n_valid, m.end - m.start + 1);
    fprintf('  - smoothed profile across window: range %.3f to %.3f\n', ...
        smooth_min, smooth_max);
    fprintf('  - raw SHAPE mean   : %.3f (SD %.3f)\n', target_mean, target_sd);
    fprintf('  - raw SHAPE median : %.3f\n', target_median);
    fprintf('  - raw SHAPE range  : %.3f to %.3f\n', target_min, target_max);
    fprintf('  - genome-wide percentile rank (median-based): %.1fth\n', pct_rank);
    fprintf('  - equivalently, top %.1f%% by accessibility\n\n', top_pct);
end

grid on;
xlabel('Nucleotide Position (nt)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('SHAPE Reactivity', 'FontSize', 12, 'FontWeight', 'bold');
title('HIV-1 Genome SHAPE Profile with Percentile Reference Lines (Top 1% to Top 90%)', 'FontSize', 13);

xlim([1 max(data.Position)]);
% Everything plotted lies on the smoothed profile, so its maximum sets the scale
ylim([0 max(smoothSHAPE, [], 'omitnan') * 1.1]);

legend(legendHandles, legendLabels, 'Location', 'northeast');

% Export at 300 dpi
exportgraphics(fig1, 'HIV1_Genome_SHAPE_Percentile_Verification_300DPI.png', 'Resolution', 300);
fprintf('Figure saved to HIV1_Genome_SHAPE_Percentile_Verification_300DPI.png (300 dpi)\n');

%% ===================== 5. DATA TIP CALLBACK =====================
% Interactive readout: clicking a point reports its position, nucleotide, raw
% reactivity and local sequence context.
dcmObj1 = datacursormode(fig1);
set(dcmObj1, 'Enable', 'on', 'UpdateFcn', {@customDataTip, data});

function txt = customDataTip(~, event_obj, dataTable)
    pos = event_obj.Position;
    x_pos = pos(1);

    [~, rowIdx] = min(abs(dataTable.Position - x_pos));

    if ~isempty(rowIdx)
        nuc = dataTable.Nucleotide(rowIdx);
        shapeVal = dataTable.SHAPE_Reactivity(rowIdx);

        startIdx = max(1, rowIdx - 5);
        endIdx = min(height(dataTable), rowIdx + 5);

        seqContext = "";
        for i = startIdx:endIdx
            if i == rowIdx
                seqContext = seqContext + "[" + dataTable.Nucleotide(i) + "]";
            else
                seqContext = seqContext + dataTable.Nucleotide(i);
            end
        end

        txt = {
            sprintf('Position: %d nt', dataTable.Position(rowIdx)), ...
            sprintf('Nucleotide: %s', nuc), ...
            sprintf('Raw SHAPE: %.3f', shapeVal), ...
            sprintf('Context (-5 to +5): %s', seqContext)
        };
    else
        txt = {sprintf('X: %d', x_pos)};
    end
end
