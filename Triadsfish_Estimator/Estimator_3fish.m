clc;
clear;
close all;

n1 = 9; % Bright
n2 = 10; % Darkness

n = n1;
dt = 1/30;

X_drift_save = zeros(n,4); % the parameters of the drift term
X_diffusion_save = zeros(n,3); % the parameters of the diffusion term

% Loading data
data_L = load('Terns_Fish_Bright_XY.mat'); % Bright: 9 datasets  
data_D = load('Terns_Fish_Dark_XY.mat'); % Darkness: 10 datasets   
data = data_L;

for i = 1:length(data.TF_B_X1)
    % Fish1
    F1_X = data.TF_B_X1{i};
    F1_Y = data.TF_B_Y1{i};

    % Fish2
    F2_X = data.TF_B_X2{i};
    F2_Y = data.TF_B_Y2{i};
    
    % Fish3
    F3_X = data.TF_B_X3{i};
    F3_Y = data.TF_B_Y3{i};

    dis12 = sqrt((F1_X - F2_X).^2 + (F1_Y - F2_Y).^2);
    dis23 = sqrt((F2_X - F3_X).^2 + (F2_Y - F3_Y).^2);
    dis13 = sqrt((F3_X - F1_X).^2 + (F3_Y - F1_Y).^2);

    d_bar_data = [dis12, dis13, dis23];
    d_bar = median(d_bar_data, 2);

    d_dot_bar =  diff(d_bar)/dt;
    d_bar = d_bar(1:end-1);
    
    interval_size = 1;


    %% Estimator
    % SSR
    Theta = constructLibrary(d_bar);

    % Kramers-Moyal Formula
    % Drift
    X_drift = zeros(size(Theta,2), 1);

    for j = 1:size(d_dot_bar,2)
        X_drift(:, :) = ssr(Theta, d_dot_bar);
    end

    X_dot_drift = Theta * X_drift;
    X_dot_drift_pred = X_dot_drift;
    X_drift_save(i,:) = X_drift';

    % Fitting
    dis_intervals = 0:interval_size:60;
    avg_dis_dot_predicted = zeros(length(dis_intervals) - 1, 1);
    avg_dis_dot_predicted2 = zeros(length(dis_intervals) - 1, 1);

    for j = 1:length(dis_intervals) - 1
        in_interval1 = (d_bar >= dis_intervals(j)) & (d_bar < dis_intervals(j + 1));
        if any(in_interval1)
            avg_dis_dot_predicted(j) = mean(X_dot_drift(in_interval1));
        else
            avg_dis_dot_predicted(j) = NaN;
        end
    end

    cell_array = {avg_dis_dot_predicted'};

    range = 1:interval_size:60;
    [beta, delta1, delta2, ci, DriftPoly, PotentialPoly] = Fitting3_1(cell_array, range, d_bar);

    % Save parameters
    betas_12(i, :) = beta;
    delta1s_12(i, :) = round(delta1, 4);
    delta2s_12(i, :) = round(delta2, 4);
    cis_12(i, :) = round(ci, 4);
    error1(i, :) = delta2 - delta1;


    % Diffusion
    Theta2 = constructLibrary2(d_bar);
    X_dot2 =  sqrt((d_dot_bar.^2) / dt^1);
    Xi_diffusion = zeros(size(Theta2,2), 1);
    for j = 1:size(X_dot2,2)
        format long;
        Xi_diffusion(:, j) = ssr2(Theta2, X_dot2(:, j));
    end
    X_dot_diff = (Theta2 * Xi_diffusion);
    X_dot_diff_pred = X_dot_diff;
    X_diffusion_save(i,:) = Xi_diffusion';

    X_dot = d_dot_bar;

    % Fitting
    for j = 1:length(dis_intervals) - 1
        in_interval1 = (d_bar >= dis_intervals(j) & d_bar < dis_intervals(j + 1));
        if any(in_interval1)
            avg_dis_dot_predicted2(j) = mean(X_dot_diff(in_interval1));
        else
            avg_dis_dot_predicted2(j) = NaN;
        end
    end

    cell_array2 = {avg_dis_dot_predicted2'};
    range = 1:interval_size:60;
    [beta, delta1, delta2, ci, Coeff1, Coeff2, Coeff3] = Fitting3_2(i, cell_array, cell_array2, range, d_bar);
    Coeff1 = round(Coeff1, 4);
    Coeff2 = round(Coeff2, 4);
    Coeff3 = round(Coeff3, 4);
    sigma_row = [Coeff1, Coeff2, Coeff3];
    sigma_row(isempty(sigma_row) | sigma_row == 0) = 0;
    sigma_row = real(sigma_row);
    sigma_12(i, :) = sigma_row;


    %% Simulation
    num_trajectories1 = 1;
    for traj = 1:num_trajectories1
        Z = zeros(1, 2*length(d_dot_bar)+1);
        Z(1) = d_bar(1);

        for ii = 1:2*length(d_dot_bar)
            if Z(ii)<0
                Z(ii) = 3;
            end

            if Z(ii)>60
                Z(ii) = 55;
            end

            % Drift term
            drift(ii) = beta * (4*(Z(ii)).^3 - 6*(delta1 + delta2)*(Z(ii)).^2 + 2*(delta1^2 + 4*delta1*delta2 + delta2^2)*(Z(ii)) - 2*delta1*delta2*(delta1 + delta2)) + ci;

            % Diffusion term
            diffusion(ii) = Coeff1 + Coeff2 * Z(ii) + Coeff3 * Z(ii).^2;

            dW = sqrt(dt) * randn(1);
            Z(ii+1) = Z(ii) + 12*drift(ii)*dt + 0.20*diffusion(ii)*dW;
        end

        Z = Z';
        Z = Z(end-18000 : end);
        Z_dot = diff(Z)/dt;
        Z = Z(1:end-1);
    end

    X_simu(:,i) = Z;
end


%% Print
print_drift = X_drift_save
print_diffusion = X_diffusion_save

% Drift
formatted_matrix = strings(size(print_drift) + [1, 0]);
for i = 1:size(print_drift, 1)
    for j = 1:size(print_drift, 2)
        val = print_drift(i, j);
        if val ~= 0 && abs(val) < 1e-3
            formatted_matrix(i, j) = sprintf('%.1e', val);
        else
            formatted_matrix(i, j) = sprintf('%06.3f', val);
        end
    end
end

% Calculate column averages for Drift
column_averages = mean(print_drift, 1);
for j = 1:size(print_drift, 2)
    if abs(column_averages(j)) < 1e-3
        formatted_matrix(end, j) = sprintf('%.1e', column_averages(j));
    else
        formatted_matrix(end, j) = sprintf('%06.3f', column_averages(j));
    end
end

% Diffusion
formatted_matrix2 = strings(size(print_diffusion) + [1, 0]);
for i = 1:size(print_diffusion, 1)
    for j = 1:size(print_diffusion, 2)
        val = print_diffusion(i, j);
        if val ~= 0 && abs(val) < 1e-3
            formatted_matrix2(i, j) = sprintf('%.1e', val);
        else
            formatted_matrix2(i, j) = sprintf('%06.3f', val);
        end
    end
end

% Calculate column averages for Diffusion
column_averages2 = mean(print_diffusion, 1);
for j = 1:size(print_diffusion, 2)
    if abs(column_averages2(j)) < 1e-3
        formatted_matrix2(end, j) = sprintf('%.1e', column_averages2(j));
    else
        formatted_matrix2(end, j) = sprintf('%06.3f', column_averages2(j));
    end
end


function [Theta, Theta_terms] = constructLibrary(X)
    syms dis;
    vars = dis;

    [n, numVars] = size(X);

    library = [];
    terms = {};

    for a1 = 0:3
        newTerm = X.^a1;
        library = [library, newTerm];

        terms{end+1} = vars^a1;
    end

    Theta = library;
    Theta_terms = terms;
end


function [Theta, Theta_terms] = constructLibrary2(X)
    syms dis;
    vars = dis;

    [n, numVars] = size(X);

    library = [];
    terms = {};

    for a1 = 0:2
        newTerm = X.^a1;
        library = [library, newTerm];

        terms{end+1} = vars^a1;
    end

    Theta = library;
    Theta_terms = terms;
end


function Xi = ssr(Theta, Y)
    alpha = 1e-4;
    max_iter = 2e+4;
    tol = 1e-8;
    Xi = lasso_init(Theta, Y, alpha, max_iter, tol);

    % Initialization
    cross_validation_scores = zeros(max_iter, 1);
    active_set = true(size(Xi)); 
    prev_Xi = Xi; 

    for iter = 1:max_iter
        mean_value = mean(abs(Xi(active_set)));
        small_values = find(abs(Xi(active_set)) < 1e-5);
        if isempty(small_values)
            break;
        else
            active_idx = find(active_set);
            min_idx_global = active_idx(small_values);

            active_set(min_idx_global) = false;
        end

        Theta_reduced = Theta(:, active_set);
        Xi_reduced = lasso_init(Theta_reduced, Y, alpha, max_iter, tol);

        Xi(active_set) = Xi_reduced;
        Xi(~active_set) = 0;

        prev_Xi = Xi;

        cross_validation_scores(iter) = cross_validation(Theta_reduced, Y);

        if iter > 1 && (cross_validation_scores(iter) / cross_validation_scores(iter - 1)) >= 1
            break;
        end

        if sum(active_set) == 0
            break;
        end
    end
end


function Xi = ssr2(Theta, Y)
    alpha = 1e-4;
    max_iter = 2e+4;
    tol = 1e-8;
    Xi = lasso_init(Theta, Y, alpha, max_iter, tol);

    % Initialization
    cross_validation_scores = zeros(max_iter, 1);
    active_set = true(size(Xi)); 
    prev_Xi = Xi; 

    for iter = 1:max_iter
        mean_value = mean(abs(Xi(active_set)));
        small_values = find(abs(Xi(active_set)) < 1e-5);
        if isempty(small_values)
            break;
        else
            active_idx = find(active_set);
            min_idx_global = active_idx(small_values);

            active_set(min_idx_global) = false;
        end

        Theta_reduced = Theta(:, active_set);
        Xi_reduced = lasso_init(Theta_reduced, Y, alpha, max_iter, tol);

        Xi(active_set) = Xi_reduced;
        Xi(~active_set) = 0;

        prev_Xi = Xi;

        cross_validation_scores(iter) = cross_validation(Theta_reduced, Y);

        if iter > 1 && (cross_validation_scores(iter) / cross_validation_scores(iter - 1)) >= 1
            break;
        end

        if sum(active_set) == 0
            break;
        end
    end
end


function Xi = lasso_init(Theta, Y, alpha, max_iter, tol)
    [n, d] = size(Theta);
    Xi = zeros(d, 1);

    Theta_norm = sum(Theta.^2, 1);

    for iter = 1:max_iter
        Xi_old = Xi;

        for j = 1:d
            r_j = Y - Theta*(Xi) + Theta(:, j)*Xi(j);
            rho = Theta(:, j)'*r_j;
            z = Theta_norm(j);

            if abs(rho)/z > alpha/z
                Xi(j) = (sign(rho) * max(0, (abs(rho) - alpha)))/z;
            else
                Xi(j) = 0;
            end
        end

        if norm(Xi - Xi_old, 2) < tol
            break;
        end
    end
end


function score = cross_validation(Theta, Y)
    K = 5;
    indices = crossvalind('Kfold', size(Y, 1), K);
    mse = zeros(K, 1);

    for i = 1:K
        test_idx = (indices == i);
        train_idx = ~test_idx;

        Theta_train = Theta(train_idx, :);
        Y_train = Y(train_idx);
        Theta_test = Theta(test_idx, :);
        Y_test = Y(test_idx);

        alpha = 1e-3;
        max_iter = 2e+4;
        tol = 1e-6;
        Xi_train = lasso_init(Theta_train, Y_train, alpha, max_iter, tol);

        Y_pred = Theta_test * Xi_train;
        mse(i) = mean((Y_test - Y_pred).^2);
    end

    score = mean(mse);
end


function result = Replace_nan(inputArray)
    % index
    nan_idx = find(isnan(inputArray));
    
    for k = 1:length(nan_idx)
        idx = nan_idx(k);

        % range
        startIdx = max(1, idx-2);
        endIdx = min(length(inputArray), idx+2);

        % delete NaN
        surroundingValues = inputArray(startIdx:endIdx);
        validValues = surroundingValues(~isnan(surroundingValues));

        % replace NaN as means
        if ~isempty(validValues)
            inputArray(idx) = mean(validValues);
        end
    end

    result = inputArray;
end


function data_no_outliers = remove_outliers(data)
    if isempty(data)
        data_no_outliers = data;
        return;
    end
    
    % Calculate the quartiles
    Q1 = quantile(data, 0.25);
    Q3 = quantile(data, 0.75);
    
    % Compute the interquartile range (IQR)
    IQR = Q3 - Q1;
    
    % Determine outliers based on IQR
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    
    % Filter out the outliers
    data_no_outliers = data(data >= lower_bound & data <= upper_bound);
end