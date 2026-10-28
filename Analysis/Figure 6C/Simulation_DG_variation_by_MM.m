%% Mismatch vs. hybridization free energy (Delta G) for 70-nt DNA duplexes
% Generates main-text Figure 6C and the mismatch-cutoff probabilities quoted
% alongside it.
%
% WHAT IT DOES
%   For each mismatch count k = 1..n, the script repeats m independent trials.
%   Each trial draws a NEW random 70-nt DNA probe with GC content in
%   [40%, 65%], builds the perfect-match target as its reverse complement,
%   introduces k mismatches at unique positions in the target, and computes
%   the probe-target hybridization free energy with NUPACK at 65 C and
%   0.5 M Na+. The distribution of Delta G at each k is then summarized as a
%   min-max band, a mean +/- SD band, and a mean line.
%   Finally, the fraction of duplexes falling below the -40 kcal/mol capture
%   cutoff is reported at each k, both empirically and under a fitted normal
%   distribution.
%
% REQUIREMENTS
%   MATLAB R2020b or later (uses exportgraphics).
%   NUPACK 4 with its MATLAB interface, plus the helper function nupackAna,
%   which must be on the MATLAB path. nupackAna wraps the NUPACK complex
%   analysis call and returns a numeric array whose element (4,1) is the
%   free energy of the probe-target complex in kcal/mol.
%   NUPACK requires a free academic licence: https://www.nupack.org
%
% HOW TO RUN
%   Place nupackAna.m on the MATLAB path, open this file and press Run, or
%   call it by name from the command window. No input arguments are needed.
%   All settings live in the "User parameters" block below.
%   Runtime scales as m * n NUPACK calls. The published settings
%   (m = 5000, n = 10) take several hours; set m = 100 for a quick test,
%   which reproduces the shape of the curve but not the published precision.
%   Outputs: one PNG at 300 dpi (name set in out_png_fig1) and a probability
%   table printed to the command window.
%
% REPRODUCIBILITY
%   The random seed is fixed below so that repeated runs give identical
%   sequences, mismatch positions and therefore identical statistics.

clear; clc;
close all;

rng(0);  % fixed seed for reproducibility

%% ---------------------- User parameters ----------------------
L            = 70;        % probe length (nt)
m            = 5000;      % trials per mismatch count (lower to ~100 for a quick test)
n            = 10;        % maximum number of mismatches tested (k = 1..n)

% Capture cutoff used in the probe design algorithm (kcal/mol)
threshold_dG = -40;

% GC content constraints for the randomly generated probe
gc_low       = 0.40;      % 40%
gc_high      = 0.65;      % 65%

% NUPACK conditions. Matched to the probe design algorithm.
% Note that the complex free energy returned here does not depend on strand
% concentration; the values are set only to keep conditions identical.
Acon         = 4e-8;      % probe strand concentration (M)
Bcon         = 1.7e-17;   % target strand concentration (M)
sodiumC      = 0.5;       % [Na+] (M)
eval_temp    = 65;        % evaluation temperature (C)

% Output figure path
out_png_fig1 = '10mm_5000_65C_.png';

%% ---------------------- Storage containers -------------------
% Summary statistics, one entry per mismatch count
mean_dG = nan(n, 1);
std_dG  = nan(n, 1);
min_dG  = nan(n, 1);
max_dG  = nan(n, 1);

% Full set of Delta G values, kept so that the empirical cutoff fraction and
% the normality check below can be computed after the loop
all_dG  = nan(m, n);

%% ---------------------- Main simulation loop -----------------
fprintf('Starting calculations at %d C (m = %d trials per mismatch count)...\n', eval_temp, m);
for k = 1:n
    dG_vals = nan(m, 1);

    fprintf('Currently calculating for mismatch count k = %d...\n', k);
    for t = 1:m
        % --- New probe for every trial, GC content in [gc_low, gc_high] ---
        probe      = rand_dna_gc(L, gc_low, gc_high);
        target_pm  = revcomp_dna(probe);

        % --- Introduce k mismatches in the target at unique positions ---
        target_mut = introduce_mismatches(target_pm, k);

        % --- NUPACK complex analysis; element (4,1) is the complex free energy ---
        results_nupack = nupackAna('dna', probe, Acon, target_mut, Bcon, eval_temp, sodiumC);
        dG_vals(t) = results_nupack(4,1);
    end

    % --- Collect statistics for this mismatch count ---
    all_dG(:, k) = dG_vals;
    mean_dG(k)   = mean(dG_vals, 'omitnan');
    std_dG(k)    = std(dG_vals, 0, 'omitnan');
    min_dG(k)    = min(dG_vals);
    max_dG(k)    = max(dG_vals);
end
fprintf('\nAll simulations complete. Generating figure.\n');

%% ---------------------- Visualization: Figure 6C ------------------------
x = (1:n)';
figure('Color','w', 'Name', 'Figure 6C: mismatch vs Delta G'); hold on; grid on;

% Polygon coordinates for the shaded bands
x_fill = [x(:)', fliplr(x(:)')];
y_fill_minmax = [min_dG(:)', fliplr(max_dG(:)')];
y_fill_sd     = [ (mean_dG(:)'-std_dG(:)') , ...
                  fliplr(mean_dG(:)'+std_dG(:)') ];

% Min-max band (blue)
hMinMax = fill(x_fill, y_fill_minmax, [0.70 0.85 1.00], ...
               'EdgeColor','none', 'FaceAlpha',0.35);

% Mean +/- SD band (orange)
hSD = fill(x_fill, y_fill_sd, [1.00 0.80 0.60], ...
           'EdgeColor','none', 'FaceAlpha',0.40);

% Mean line (black)
hMean = plot(x, mean_dG, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'Color', 'k');

xlabel('Number of mismatches', 'FontSize', 16);
ylabel('\DeltaG (kcal/mol)', 'FontSize', 16);
legend([hMinMax, hSD, hMean], {'Min-max range','Mean \pm SD',['Mean (', num2str(eval_temp), '^\circC)']}, ...
       'Location','best', 'Box','off', 'FontSize', 16);
xlim([1 n]);
set(gca, 'FontSize', 16);

% Export figure
f1 = gcf;
exportgraphics(f1, out_png_fig1, 'Resolution', 300);
fprintf('Figure saved to %s (300 dpi)\n', out_png_fig1);

%% ---------------------- Fraction below the capture cutoff -------------
% Two estimates are reported at each mismatch count:
%   Empirical  - the observed fraction of trials with Delta G below the cutoff.
%                This makes no distributional assumption and is the primary value.
%   Normal fit - the normal CDF evaluated at the cutoff using the sample mean
%                and SD. Reported for comparison; close agreement with the
%                empirical fraction indicates the normal approximation is
%                adequate in the tail that matters here.
fprintf('\n--- Fraction of duplexes with Delta G < %.1f kcal/mol at %d C ---\n', threshold_dG, eval_temp);
fprintf('%4s  %12s  %12s  %10s\n', 'k', 'Empirical', 'Normal fit', 'Difference');
emp_frac  = nan(n, 1);
norm_frac = nan(n, 1);
for k = 1:n
    vals = all_dG(~isnan(all_dG(:,k)), k);
    emp_frac(k)  = sum(vals < threshold_dG) / numel(vals);
    norm_frac(k) = normcdf(threshold_dG, mean_dG(k), std_dG(k));
    fprintf('%4d  %11.2f%%  %11.2f%%  %9.2f%%\n', ...
        k, emp_frac(k)*100, norm_frac(k)*100, (emp_frac(k)-norm_frac(k))*100);
end
fprintf('-----------------------------------------------------------\n\n');

%% ---------------------- Normality check -------------------------------
% Tests whether the Delta G values at each mismatch count are consistent with
% a normal distribution, which is the assumption behind the normal-fit column
% above. The Lilliefors test is used because the mean and SD are estimated
% from the data rather than specified in advance. With m in the thousands the
% test is sensitive to small departures, so the skewness and excess kurtosis
% are printed alongside it as effect sizes: values near zero indicate that any
% departure is slight. Requires the Statistics and Machine Learning Toolbox;
% the check is skipped if lillietest is unavailable.
if exist('lillietest', 'file') == 2
    fprintf('--- Normality check on the Delta G distribution at each k ---\n');
    fprintf('%4s  %10s  %10s  %12s  %10s\n', 'k', 'Skewness', 'Ex.Kurt', 'Lilliefors p', 'Normal?');
    for k = 1:n
        vals = all_dG(~isnan(all_dG(:,k)), k);
        sk = skewness(vals);
        ek = kurtosis(vals) - 3;
        % lillietest caps its tabulated p-values, so p is returned as a bound
        [h, p] = lillietest(vals, 'Alpha', 0.05);
        if h == 0, verdict = 'yes'; else, verdict = 'no'; end
        fprintf('%4d  %10.3f  %10.3f  %12.4f  %10s\n', k, sk, ek, p, verdict);
    end
    fprintf('-----------------------------------------------------------\n\n');
else
    fprintf('lillietest not available; skipping the normality check.\n');
    fprintf('The empirical fractions above do not depend on this check.\n\n');
end

%% ---------------------- Helper functions ---------------------
% NOTE: these functions must appear only once, at the end of the script.

function seq = rand_dna_gc(L, gc_low, gc_high)
    % Random DNA sequence of length L with GC content drawn uniformly from
    % [gc_low, gc_high].
    if ~(gc_low >= 0 && gc_high <= 1 && gc_low <= gc_high)
        error('gc_low/high must satisfy 0 <= gc_low <= gc_high <= 1.');
    end
    f_gc     = gc_low + (gc_high - gc_low) * rand();
    gc_count = round(L * f_gc);

    pos_gc = false(1, L);
    if gc_count > 0
        pos_gc(randperm(L, gc_count)) = true;
    end
    pos_at = ~pos_gc;

    seq = repmat('A', 1, L);
    gc_chars = 'GC';
    idx_gc   = find(pos_gc);
    if ~isempty(idx_gc)
        seq(idx_gc) = gc_chars(randi(2, 1, numel(idx_gc)));
    end

    at_chars = 'AT';
    idx_at   = find(pos_at);
    if ~isempty(idx_at)
        seq(idx_at) = at_chars(randi(2, 1, numel(idx_at)));
    end
end

function rc = revcomp_dna(seq)
    % Reverse complement of a DNA sequence.
    seq = upper(seq);
    cmap = containers.Map({'A','C','G','T'}, {'T','G','C','A'});
    comp = arrayfun(@(c) cmap(c), seq);
    rc   = fliplr(comp);
end

function out = introduce_mismatches(seq, k)
    % Replaces the base at k unique random positions with one of the three
    % other bases, so every substitution is a genuine mismatch.
    seq = upper(seq);
    L   = numel(seq);

    if k < 0 || k > L
        error('k must be between 0 and sequence length.');
    end
    if k == 0
        out = seq;
        return;
    end

    pos = randperm(L, k);
    bases = 'ACGT';
    out = seq;
    for i = 1:k
        p   = pos(i);
        cur = out(p);
        cand = bases(bases ~= cur);
        out(p) = cand(randi(3));
    end
end
