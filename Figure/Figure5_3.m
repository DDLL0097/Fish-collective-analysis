clc
clear all
close all

% loead data
for i = 1:9
    % Bright condition
    fileName = sprintf('Est_diff_L_%d.mat', i);
    filePath = fullfile('Triads_Bright', fileName);
    DatA = load(filePath);
    fileName2 = sprintf('Est_drift_L_%d.mat', i);
    filePath2 = fullfile('Triads_Bright', fileName2);
    DatB = load(filePath2);
    Data_B_Est_Difusion{i} = DatA.X_dot_diff;
    Data_B_Est_Drift{i} = DatB.X_dot_drift;
end

for i = 1:10
    % Dark condition
    fileName3 = sprintf('Est_diff_D_%d.mat', i);
    filePath3 = fullfile('Triads_Dark', fileName3);
    DatC = load(filePath3);
    fileName4 = sprintf('Est_drift_D_%d.mat', i);
    filePath4 = fullfile('Triads_Dark', fileName4);
    DatD = load(filePath4);
    Data_D_Est_Difusion{i} = DatC.X_dot_diff;
    Data_D_Est_Drift{i} = DatD.X_dot_drift;
end

% position
TF_D = load('Triads_Fish_Dark.mat').TF_D;
TF_B = load('Triads_Fish_Bright.mat').TF_B;

% Filtering
order = 4;
wc = 0.1;
[b, a] = butter(order, wc, 'low');


%% ========================================================
% Bright
for i=1:length(TF_B)
    X1 = TF_B{i}.X1; Y1 = TF_B{i}.Y1;
    X2 = TF_B{i}.X2; Y2 = TF_B{i}.Y2;
    X3 = TF_B{i}.X2; Y3 = TF_B{i}.Y2;

    X1 = fillmissing(X1, 'next'); Y1 = fillmissing(Y1, 'next');
    X2 = fillmissing(X2, 'next'); Y2 = fillmissing(Y2, 'next');
    X3 = fillmissing(X3, 'next'); Y3 = fillmissing(Y3, 'next');

    X1f = filtfilt(b, a, X1); Y1f = filtfilt(b, a, Y1);
    X2f = filtfilt(b, a, X2); Y2f = filtfilt(b, a, Y2);
    X3f = filtfilt(b, a, X3); Y3f = filtfilt(b, a, Y3);

    X1 = X1f; Y1 = Y1f;
    X2 = X2f; Y2 = Y2f;
    X3 = X3f; Y3 = Y3f;

    D1 = sqrt((X1 - X2).^2+(Y1 - Y2).^2);
    D2 = sqrt((X2 - X3).^2+(Y2 - Y3).^2);
    D3 = sqrt((X1 - X3).^2+(Y1 - Y3).^2);

    d_bar_data = [D1, D2, D3];
    d_bar = median(d_bar_data, 2);

    Data_B_D{i} = d_bar(1:end-1);
end

% Dark
for i=1:length(TF_D)
    X1 = TF_D{i}.X1; Y1 = TF_D{i}.Y1;
    X2 = TF_D{i}.X2; Y2 = TF_D{i}.Y2;
    X3 = TF_D{i}.X2; Y3 = TF_D{i}.Y2;

    X1 = fillmissing(X1, 'next'); Y1 = fillmissing(Y1, 'next');
    X2 = fillmissing(X2, 'next'); Y2 = fillmissing(Y2, 'next');
    X3 = fillmissing(X3, 'next'); Y3 = fillmissing(Y3, 'next');

    X1f = filtfilt(b, a, X1); Y1f = filtfilt(b, a, Y1);
    X2f = filtfilt(b, a, X2); Y2f = filtfilt(b, a, Y2);
    X3f = filtfilt(b, a, X3); Y3f = filtfilt(b, a, Y3);

    X1 = X1f; Y1 = Y1f;
    X2 = X2f; Y2 = Y2f;
    X3 = X3f; Y3 = Y3f;

    D1 = sqrt((X1 - X2).^2+(Y1 - Y2).^2);
    D2 = sqrt((X2 - X3).^2+(Y2 - Y3).^2);
    D3 = sqrt((X1 - X3).^2+(Y1 - Y3).^2);

    d_bar_data = [D1, D2, D3];
    d_bar = median(d_bar_data, 2);

    Data_D_D{i} = d_bar(1:end-1);
end


%======= plot of average function =========================
%======= BRIGHT CONDITION =================================
% Initialize
integrals_Y_squared = zeros(1, length(Data_B_Est_Difusion));
integrals_X_squared = zeros(1, length(Data_B_Est_Difusion));

for r = 2:length(Data_B_Est_Difusion)
    X = Data_B_Est_Difusion{r};
    Y = Data_B_Est_Drift{r};
    D = Data_B_D{r};
    
    % Ensure data alignment and remove NaNs
    valid_idx = ~isnan(Y) & ~isnan(X) & ~isnan(D);
    D = D(valid_idx);
    Y = Y(valid_idx);
    X = X(valid_idx);
    
    % Sort D, and rearrange Y and X accordingly
    [D, sort_idx] = sort(D);
    Y = Y(sort_idx);
    X = X(sort_idx);
    
    % Compute the integral of Y^2 and X^2 over D
    integrals_Y_squared(r) = trapz(D, Y.^2);
    integrals_X_squared(r) = trapz(D, X.^2);
end

% Display the integrals
disp('Integrals of Y^2 over D for each dataset:');
disp(integrals_Y_squared);

disp('Integrals of X^2 over D for each dataset:');
disp(integrals_X_squared);

% Parameters for binning
epsilon = 0.1; % Bin width
bin_edges = min(cellfun(@min, Data_B_D)) : epsilon : max(cellfun(@max, Data_B_D)); % Define bin edges
bin_centers = bin_edges(1:end-1) + epsilon / 2; % Compute bin centers


%% Process for Y
% Initialize
all_Y = [];
all_D_Y = [];

% Collect all D and Y values
for r = 1:length(Data_B_Est_Drift)
    all_Y = [all_Y; Data_B_Est_Drift{r}];
    all_D_Y = [all_D_Y; Data_B_D{r}];
end

% Group Y values into bins based on D
bin_indices_Y = discretize(all_D_Y, bin_edges); % Assign each D to a bin

mean_Y_bins = accumarray(bin_indices_Y(~isnan(bin_indices_Y)), all_Y(~isnan(bin_indices_Y)), ...
    [length(bin_centers), 1], @mean, NaN); % Mean Y for each bin
std_Y_bins = accumarray(bin_indices_Y(~isnan(bin_indices_Y)), all_Y(~isnan(bin_indices_Y)), ...
    [length(bin_centers), 1], @std, NaN); % Std Y for each bin

% Plot the averaged data for Y with shaded regions
figure;

% Shaded region for std (Y)
subplot(1,2,1);
fill([bin_centers, fliplr(bin_centers)], ...
    [mean_Y_bins + std_Y_bins; flipud(mean_Y_bins - std_Y_bins)], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter blue with transparency

hold on;

% Plot mean line (Y)
plot(bin_centers, mean_Y_bins, '-b', 'LineWidth', 2); % Darker blue for mean line

xlabel('D');
ylabel('Average Y');
title('Binned Average of Y vs. D');


%% Process for X
% Initialize
all_X = [];
all_D_X = [];

% Collect all D and X values
for r = 1:length(Data_B_Est_Difusion)
    all_X = [all_X; Data_B_Est_Difusion{r}];
    all_D_X = [all_D_X; Data_B_D{r}];
end

% Group X values into bins based on D
bin_indices_X = discretize(all_D_X, bin_edges); % Assign each D to a bin

mean_X_bins = accumarray(bin_indices_X(~isnan(bin_indices_X)), all_X(~isnan(bin_indices_X)), ...
    [length(bin_centers), 1], @mean, NaN); % Mean X for each bin
std_X_bins = accumarray(bin_indices_X(~isnan(bin_indices_X)), all_X(~isnan(bin_indices_X)), ...
    [length(bin_centers), 1], @std, NaN); % Std X for each bin

% Plot the averaged data for X with shaded regions
subplot(1,2,2);
fill([bin_centers, fliplr(bin_centers)], ...
    [mean_X_bins + std_X_bins; flipud(mean_X_bins - std_X_bins)], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter blue with transparency

hold on;

% Plot mean line (X)
plot(bin_centers, mean_X_bins, '-b', 'LineWidth', 2); % Darker blue for mean line

xlabel('D');
ylabel('Average X');
title('Binned Average of X vs. D');

all_YBright = all_Y;
all_XBright = all_X;


%======= plot of average function =========================
%======= DARK CONDITION ===================================
% Initialize
integrals_Y_squared_D = zeros(1, length(Data_D_Est_Difusion));
integrals_X_squared_D = zeros(1, length(Data_D_Est_Difusion));

for r = 1:length(Data_D_Est_Difusion)
    X = Data_D_Est_Difusion{r};
    Y = Data_D_Est_Drift{r};
    D = Data_D_D{r};
    
    % Ensure data alignment and remove NaNs
    valid_idx = ~isnan(Y) & ~isnan(X) & ~isnan(D);
    D = D(valid_idx);
    Y = Y(valid_idx);
    X = X(valid_idx);
    
    % Sort D, and rearrange Y and X accordingly
    [D, sort_idx] = sort(D);
    Y = Y(sort_idx);
    X = X(sort_idx);
    
    % Compute the integral of Y^2 and X^2 over D
    integrals_Y_squared_D(r) = trapz(D, Y.^2);
    integrals_X_squared_D(r) = trapz(D, X.^2);
end

% Display the integrals
disp('Integrals of Y^2 over D for each dataset:');
disp(integrals_Y_squared);

disp('Integrals of X^2 over D for each dataset:');
disp(integrals_X_squared);

% Parameters for binning
epsilon = 0.1; % Bin width
bin_edges = min(cellfun(@min, Data_D_D)) : epsilon : max(cellfun(@max, Data_D_D)); % Define bin edges
bin_centers2 = bin_edges(1:end-1) + epsilon / 2; % Compute bin centers


%% Process for Y
% Initialize
all_Y = [];
all_D_Y = [];

% Collect all D and Y values
for r = 1:length(Data_D_Est_Drift)
    all_Y = [all_Y; Data_D_Est_Drift{r}];
    all_D_Y = [all_D_Y; Data_D_D{r}];
end

% Group Y values into bins based on D
bin_indices_Y = discretize(all_D_Y, bin_edges); % Assign each D to a bin

mean_Y_bins2 = accumarray(bin_indices_Y(~isnan(bin_indices_Y)), all_Y(~isnan(bin_indices_Y)), ...
    [length(bin_centers2), 1], @mean, NaN); % Mean Y for each bin
std_Y_bins2 = accumarray(bin_indices_Y(~isnan(bin_indices_Y)), all_Y(~isnan(bin_indices_Y)), ...
    [length(bin_centers2), 1], @std, NaN); % Std Y for each bin

% Plot the averaged data for Y with shaded regions
figure;

% Shaded region for std (Y)
subplot(1,2,1);
fill([bin_centers2, fliplr(bin_centers2)], ...
    [mean_Y_bins2 + std_Y_bins2; flipud(mean_Y_bins2 - std_Y_bins2)], ...
    [1, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter blue with transparency

hold on;

% Plot mean line (Y)
plot(bin_centers2, mean_Y_bins2, '-r', 'LineWidth', 2); % Darker blue for mean line

xlabel('D');
ylabel('Average Y');
title('Binned Average of Y vs. D');


%% Process for X
% Initialize
all_X = [];
all_D_X = [];

% Collect all D and X values
for r = 1:length(Data_D_Est_Difusion)
    all_X = [all_X; Data_D_Est_Difusion{r}];
    all_D_X = [all_D_X; Data_D_D{r}];
end

% Group X values into bins based on D
bin_indices_X = discretize(all_D_X, bin_edges); % Assign each D to a bin

mean_X_bins2 = accumarray(bin_indices_X(~isnan(bin_indices_X)), all_X(~isnan(bin_indices_X)), ...
    [length(bin_centers2), 1], @mean, NaN); % Mean X for each bin
std_X_bins2 = accumarray(bin_indices_X(~isnan(bin_indices_X)), all_X(~isnan(bin_indices_X)), ...
    [length(bin_centers2), 1], @std, NaN); % Std X for each bin

% Plot the averaged data for X with shaded regions
subplot(1,2,2);
fill([bin_centers2, fliplr(bin_centers2)], ...
    [mean_X_bins2 + std_X_bins2; flipud(mean_X_bins2 - std_X_bins2)], ...
    [1, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter blue with transparency

hold on;

% Plot mean line (X)
plot(bin_centers2, mean_X_bins2, '-r', 'LineWidth', 2); % Darker blue for mean line

xlabel('D');
ylabel('Average X');
title('Binned Average of X vs. D');

all_YDark = all_Y; % drift 
all_XDark = all_X; % diffusion


%=================== Final Figure =====================================
fsize = 22;

% First Figure: Y-Axis Data
figure;
set(gcf, 'Color', 'w', 'Position', [100, 100, 603, 471]); % Set figure background to white and figure size

fill([bin_centers, fliplr(bin_centers)], ...
    [mean_Y_bins + std_Y_bins; flipud(mean_Y_bins - std_Y_bins)], ...
    [1, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter blue with transparency
hold on;
h1 = plot(bin_centers, mean_Y_bins, '-r', 'LineWidth', 2); % Darker blue for mean line

fill([bin_centers2, fliplr(bin_centers2)], ...
    [mean_Y_bins2 + std_Y_bins2; flipud(mean_Y_bins2 - std_Y_bins2)], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter red with transparency
h2 = plot(bin_centers2, mean_Y_bins2, '-b', 'LineWidth', 2); % Darker red for mean line
grid off;
set(gca, 'GridLineStyle', '--'); % Dashed grid lines
set(gca, 'Color', 'w'); % Set axis background to white
set(gca, 'FontName', 'Helvetica', 'FontSize', fsize);
axis([0 60 -0.8 0.8])
legend([h1, h2], {'Bright', 'Dark'}, ...
       'FontName', 'Helvetica', 'FontSize', fsize, 'Location', 'best', 'Box', 'off'); 

% Second Figure: X-Axis Data
figure;
set(gcf, 'Color', 'w', 'Position', [100, 100, 603, 471]); % Set figure background to white and figure size

fill([bin_centers, fliplr(bin_centers)], ...
    [mean_X_bins + std_X_bins; flipud(mean_X_bins - std_X_bins)], ...
    [1, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter blue with transparency
hold on;
h1 = plot(bin_centers, mean_X_bins, '-r', 'LineWidth', 2); % Darker blue for mean line

fill([bin_centers2, fliplr(bin_centers2)], ...
    [mean_X_bins2 + std_X_bins2; flipud(mean_X_bins2 - std_X_bins2)], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Lighter red with transparency
h2 = plot(bin_centers2, mean_X_bins2, '-b', 'LineWidth', 2); % Darker red for mean line

grid off;
set(gca, 'GridLineStyle', '--'); % Dashed grid lines
set(gca, 'Color', 'w'); % Set axis background to white
set(gca, 'FontName', 'Helvetica', 'FontSize', fsize);
axis([0 60 0 80])
legend([h1, h2], {'Bright', 'Dark'}, ...
       'FontName', 'Helvetica', 'FontSize', fsize, 'Location', 'best', 'Box', 'off'); 


%=================== Integral of the potential Figure =================
integral_Y = cumtrapz(bin_centers, mean_Y_bins);
integral_Y2 = cumtrapz(bin_centers2, mean_Y_bins2);

% Plot the results
figure;
set(gcf, 'Color', 'w', 'Position', [100, 100, 603, 471]); % Set figure background to white and figure size
h1 = plot(bin_centers, -integral_Y, '-r', 'LineWidth', 2); % Integral plot
hold on
h2 = plot(bin_centers2,-integral_Y2, '-b', 'LineWidth', 2); % Darker red for mean line

grid off;
set(gca, 'GridLineStyle', '--'); % Dashed grid lines
set(gca, 'Color', 'w'); % Set axis background to white
set(gca, 'FontName', 'Helvetica', 'FontSize', fsize);
axis([0 60 -3 8])
legend([h1, h2], {'Bright', 'Dark'}, ...
       'FontName', 'Helvetica', 'FontSize', fsize, 'Location', 'best', 'Box', 'off');