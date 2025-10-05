function [beta, delta1, delta2, ci, DriftPoly, PotentialPoly] = Fitting3_1(Drift, range, x1)
    % Initialize cell arrays to collect drift and diffusion polynomial coefficients
    DriftPoly = cell(length(Drift), 1);
    PotentialPoly = cell(length(Drift), 1);
    
    % Define the function to fit for drift
    fit_func_drift = @(params, x1) 1 * (params(1) * (4*x1.^3 - 6*(params(2) + params(3))*x1.^2 + 2*(params(2)^2 + 4*params(2)*params(3) + params(3)^2)*x1 - 2*params(2)*params(3)*(params(2) + params(3))) + params(4));

    % Initialize arrays to collect drift and diffusion data
    all_drift_data = [];
    all_range_data = [];
    all_diff_data = [];
    all_range_data2 = [];
    
    % Collect and fit the data
    for i = 1:length(Drift)

        % Fit drift data
        validIdx = ~isnan(Drift{i});
        current_drift = Drift{i}(validIdx);
        all_drift_data = [all_drift_data; current_drift];
        all_range_data = [all_range_data; range(validIdx)];
    end
    
    % Initial guess
    beta_initial = max(abs(all_drift_data)) / (max(all_range_data)^3); % Scale using the cube of the data range
    ci_initial = mean(all_drift_data); % initial guess

    % local min and max
    isMin = islocalmin(all_drift_data); % Find local minima
    isMax = islocalmax(all_drift_data); % Find local maxima
    
    min_locs = all_range_data(isMin(:)); % x-coordinates of local minima
    max_locs = all_range_data(isMax(:)); % x-coordinates of local maxima
    min_vals = all_drift_data(isMin(:)); % y-coordinates of local minima
    max_vals = all_drift_data(isMax(:)); % y-coordinates of local maxima
    
    % merge
    min_locs = min_locs(:);
    max_locs = max_locs(:);
    min_vals = min_vals(:);
    max_vals = max_vals(:); 

    extrema_locs = [min_locs; max_locs]; % Concatenated x-coordinates
    extrema_vals = [min_vals; max_vals]; % Concatenated y-coordinates
    
    [extrema_locs, unique_idx] = unique(extrema_locs, 'stable');
    extrema_vals = extrema_vals(unique_idx);
    
    if isempty(extrema_locs) || isempty(extrema_vals)
        warning('isempty');
    end

    % Sort by x-coordinates to determine delta1 and delta2
    [sorted_locs, sort_idx] = sort(extrema_locs);
    sorted_vals = extrema_vals(sort_idx);

    if length(sorted_locs) < 2
        delta1_initial = sorted_locs;
        delta2_initial = sorted_locs;        
    else

        delta1_initial = sorted_locs(1); % delta1 is on the left
        delta2_initial = sorted_locs(2); % delta2 is on the right
    end

    % Perform the non-linear least squares fitting for drift with bounds
    if length(sorted_locs) < 2 || isempty(sorted_locs(:,:))
        if isempty(sorted_locs)
            sorted_locs = 60;
            delta1_initial = sorted_locs;
            delta2_initial = sorted_locs;
        else

            delta1_initial = 0;
            delta2_initial = sorted_locs;
        end

        initial_guess = [beta_initial, delta1_initial, delta2_initial, ci_initial];
        lb = [-Inf, 0, 0, -Inf];
        ub = [Inf, 60, 60, Inf];
        options = optimset('Display', 'off');
        [params_est, resnorm] = lsqcurvefit(fit_func_drift, initial_guess, all_range_data, all_drift_data, lb, ub, options);

        % Extract the estimated parameters
        beta = params_est(1);
        delta1 = params_est(2);
        delta2 = params_est(3);
        ci = params_est(4);
    else

        initial_guess = [beta_initial, delta1_initial, delta2_initial, ci_initial];
        lb = [-Inf, 0, 0, -Inf];
        ub = [Inf, 60, 60, Inf];
        options = optimset('Display', 'off');
        [params_est, resnorm] = lsqcurvefit(fit_func_drift, initial_guess, all_range_data, all_drift_data, lb, ub, options);

        % Extract the estimated parameters
        beta = params_est(1);
        delta1 = params_est(2);
        delta2 = params_est(3);
        ci = params_est(4);
    end

    % % Plot
    % colorData1 = Drift{1};
    % colors = jet(length(colorData1));
    % 
    % % Drift
    % figure;
    % hold on;
    % for i = 1:length(Drift)
    %     scatter(range', Drift{i}, 36, colorData1, 'filled');
    %     disp('1:');
    %     disp(params_est)
    %     plot(range, fit_func_drift(params_est, range), '-', 'Color', [0.9, 0.1, 0.1], 'LineWidth', 1);
    %     %plot(range, polyval(DriftPoly{i}, range), '-', 'Color', [0.9, 0.1, 0.1], 'LineWidth', 1);
    % end
    % xlabel('Euclidean Distance (cm)', 'FontSize', 20);
    % ylabel('Drift Term', 'FontSize', 20);
    % colorbar;
    % caxis([-0.3 0.6]);
    % %ylim([-20 60]);
end