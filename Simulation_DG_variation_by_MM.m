%% Mismatch vs. ΔG (min–max band + mean±SD band + mean line), font size 14
% - Each trial generates a NEW 70-nt probe with GC content in [40%, 65%]
% - Perfect-match target = reverse complement of probe
% - Evaluates at a single temperature: 65 C
% - Figure 1: 65C Details (min-max, mean±SD bands)
clear; clc;
close all;

%% ---------------------- User parameters ----------------------
L            = 70;        % probe length (nt)
m            = 5000;      % trials per mismatch count (Recommend lowering to around 100 for initial testing)
n            = 10;        % maximum number of mismatches to test (k = 1..n)

% Threshold for probability calculation (kcal/mol)
threshold_dG = -41;       

% GC content constraints for the probe
gc_low       = 0.40;      % 40%
gc_high      = 0.65;      % 65%

% NUPACK-related conditions
Acon         = 1;         % strand A concentration (M) for probe
Bcon         = 1;         % strand B concentration (M) for target
sodiumC      = 0.5;       % [Na+] (M)
eval_temp    = 65;        % Temperature for evaluation (C)

% Output figure paths
out_png_fig1 = '10mm_5000_65C_.png';

%% ---------------------- Storage containers -------------------
% Store as an (n mismatch_counts) x 1 array
mean_dG = nan(n, 1);
std_dG  = nan(n, 1);
min_dG  = nan(n, 1);
max_dG  = nan(n, 1);

%% ---------------------- Main simulation loop -----------------
fprintf('Starting calculations for temperature %d C...\n', eval_temp);
for k = 1:n
    % Temporary data array of (m trials) for the current mismatch k
    dG_vals = nan(m, 1);
    
    fprintf('Currently calculating for mismatch count k = %d...\n', k);
    for t = 1:m
        % --- New probe for every trial with GC in [gc_low, gc_high] ---
        probe      = rand_dna_gc(L, gc_low, gc_high);   
        target_pm  = revcomp_dna(probe);                
        
        % --- Introduce k mismatches in the target (unique positions) ---
        target_mut = introduce_mismatches(target_pm, k);
        
        % --- NUPACK calculation ---
        results_nupack = nupackAna('dna', probe, Acon, target_mut, Bcon, eval_temp, sodiumC);
        dG_vals(t) = results_nupack(4,1);
    end
    
    % --- Collect stats for this mismatch count k ---
    mean_dG(k) = mean(dG_vals, 'omitnan');
    std_dG(k)  = std(dG_vals, 0, 'omitnan');
    min_dG(k)  = min(dG_vals);
    max_dG(k)  = max(dG_vals);
end
fprintf('\nAll simulations complete! Generating graph.\n');

%% ---------------------- Visualization: Figure 1 (65℃) ------------------------
x = (1:n)';
figure('Color','w', 'Name', 'Figure 1: 65C Details'); hold on; grid on;

% Build polygon coordinates
x_fill = [x(:)', fliplr(x(:)')];                               
y_fill_minmax = [min_dG(:)', fliplr(max_dG(:)')];              
y_fill_sd     = [ (mean_dG(:)'-std_dG(:)') , ...               
                  fliplr(mean_dG(:)'+std_dG(:)') ];            

% Min–max band (blue)
hMinMax = fill(x_fill, y_fill_minmax, [0.70 0.85 1.00], ...
               'EdgeColor','none', 'FaceAlpha',0.35);

% Mean ± SD band (orange)
hSD = fill(x_fill, y_fill_sd, [1.00 0.80 0.60], ...
           'EdgeColor','none', 'FaceAlpha',0.40);

% Mean line (black)
hMean = plot(x, mean_dG, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'Color', 'k');

xlabel('Number of mismatches', 'FontSize', 16);
ylabel('\DeltaG (kcal/mol)', 'FontSize', 16);
legend([hMinMax, hSD, hMean], {'Min–max range','Mean \pm SD',['Mean (', num2str(eval_temp), '^\circC)']}, ...
       'Location','best', 'Box','off', 'FontSize', 16);
xlim([1 10]);
set(gca, 'FontSize', 16);

% Export Figure 1
f1 = gcf;
exportgraphics(f1, out_png_fig1, 'Resolution', 300);
fprintf('Figure 1 saved to %s (300 DPI)\n', out_png_fig1);

%% ---------------------- Probability Calculation (Dynamic Threshold) -------------
% Calculate the probability of dG < threshold_dG assuming a normal distribution
fprintf('\n--- Probability of dG < %.1f kcal/mol at %d degrees C ---\n', threshold_dG, eval_temp);
for k = 1:n
    % normcdf calculates the normal cumulative distribution function
    prob = normcdf(threshold_dG, mean_dG(k), std_dG(k));
    fprintf('Mismatch k = %2d: Probability = %.4f (%.2f%%)\n', k, prob, prob * 100);
end
fprintf('----------------------------------------------------------\n\n');

%% ---------------------- Helper functions ---------------------
% NOTE: Do not duplicate these functions! They must appear only once at the end of the script.

function seq = rand_dna_gc(L, gc_low, gc_high)
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
    seq = upper(seq);
    cmap = containers.Map({'A','C','G','T'}, {'T','G','C','A'});
    comp = arrayfun(@(c) cmap(c), seq);  
    rc   = fliplr(comp);                 
end

function out = introduce_mismatches(seq, k)
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