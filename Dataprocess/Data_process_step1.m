clc;
clear;
close all;

n1 = 9;
n2 = 10;
dt = 1/30;

% Single
for i = 1:n2
    filename = sprintf('Single_L_%d.csv', i);
    A = readmatrix(filename);
    time = A(1:18000,1);
    
    % Rescaling x and y;
    % inch to cm
    unit_c = 2.54;
    % pixel to cm
    alpha = unit_c*2400/(2048*100);                         %scaling number
    
    F1_X = A(1:18000,3);
    F1_Y = A(1:18000,4);
    
    % replace nan values
    F1_X = Replace_nan(F1_X);
    F1_Y = Replace_nan(F1_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F1_X = filtfilt(b,a,F1_X);
    F1_Y = filtfilt(b,a,F1_Y);
    
    F1_X = (alpha*F1_X) ;
    F1_Y = (alpha*F1_Y) ;
    
    minValue_F1_X = min(F1_X);
    minValue_F1_Y = min(F1_Y);
    maxValue_F1_X = max(F1_X);
    maxValue_F1_Y = max(F1_Y);
    beta1 = -((maxValue_F1_X+minValue_F1_X)/2);             %scaling number
    beta2 = -((maxValue_F1_Y+minValue_F1_Y)/2);    
    
    F1_X = (F1_X + beta1) ;
    F1_Y = (F1_Y + beta2) ;
    
    SF_B{i}.X = F1_X;
    SF_B{i}.Y = F1_Y;
end
save('Single_Fish_Bright.mat',"SF_B")

for i = 1:n2
    filename = sprintf('Single_D_%d.csv', i);
    A = readmatrix(filename);
    time = A(1:18000,1);
    
    % Rescaling x and y;
    % inch to cm
    unit_c = 2.54;
    % pixel to cm
    alpha = unit_c*2400/(2048*100);                         %scaling number
    
    F1_X = A(1:18000,3);
    F1_Y = A(1:18000,4);
    
    % replace nan values
    F1_X = Replace_nan(F1_X);
    F1_Y = Replace_nan(F1_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F1_X = filtfilt(b,a,F1_X);
    F1_Y = filtfilt(b,a,F1_Y);
    
    F1_X = (alpha*F1_X) ;
    F1_Y = (alpha*F1_Y) ;
    
    minValue_F1_X = min(F1_X);
    minValue_F1_Y = min(F1_Y);
    maxValue_F1_X = max(F1_X);
    maxValue_F1_Y = max(F1_Y);
    beta1 = -((maxValue_F1_X+minValue_F1_X)/2);             %scaling number
    beta2 = -((maxValue_F1_Y+minValue_F1_Y)/2);    
    
    F1_X = (F1_X + beta1) ;
    F1_Y = (F1_Y + beta2) ;
    
    SF_D{i}.X = F1_X;
    SF_D{i}.Y = F1_Y;
end
save('Single_Fish_Dark.mat',"SF_D")

% Pairs
for i = 1:n1
    filename = sprintf('Pair_L_%d.csv', i);
    A = readmatrix(filename);
    time = A(1:18000,1);
    
    % Rescaling x and y;
    % inch to cm
    unit_c = 2.54;
    % pixel to cm
    alpha = unit_c*2400/(2048*100);                         %scaling number
    
    F1_X = A(1:18000,3);
    F1_Y = A(1:18000,4);
    
    % replace nan values
    F1_X = Replace_nan(F1_X);
    F1_Y = Replace_nan(F1_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F1_X = filtfilt(b,a,F1_X);
    F1_Y = filtfilt(b,a,F1_Y);
    
    F1_X = (alpha*F1_X) ;
    F1_Y = (alpha*F1_Y) ;
    
    minValue_F1_X = min(F1_X);
    minValue_F1_Y = min(F1_Y);
    maxValue_F1_X = max(F1_X);
    maxValue_F1_Y = max(F1_Y);
    beta1 = -((maxValue_F1_X+minValue_F1_X)/2);             %scaling number
    beta2 = -((maxValue_F1_Y+minValue_F1_Y)/2);    
    
    F1_X = (F1_X + beta1) ;
    F1_Y = (F1_Y + beta2) ;
    
    F2_X = A(1:18000,5);
    F2_Y = A(1:18000,6);
    
    % replace nan values
    F2_X = Replace_nan(F2_X);
    F2_Y = Replace_nan(F2_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F2_X = filtfilt(b,a,F2_X);
    F2_Y = filtfilt(b,a,F2_Y);
    
    F2_X = (alpha*F2_X) ;
    F2_Y = (alpha*F2_Y) ;
    
    minValue_F2_X = min(F2_X);
    minValue_F2_Y = min(F2_Y);
    maxValue_F2_X = max(F2_X);
    maxValue_F2_Y = max(F2_Y);
    beta3 = -((maxValue_F2_X+minValue_F2_X)/2);             %scaling number
    beta4 = -((maxValue_F2_Y+minValue_F2_Y)/2); 
    
    F2_X = (F2_X + beta3) ;
    F2_Y = (F2_Y + beta4) ;
    
    PF_B{i}.X1 = F1_X;
    PF_B{i}.Y1 = F1_Y;
    PF_B{i}.X2 = F2_X;
    PF_B{i}.Y2 = F2_Y;
end
save('Pairs_Fish_Bright.mat',"PF_B")

for i = 1:n2
    filename = sprintf('Pair_D_%d.csv', i);
    A = readmatrix(filename);
    time = A(1:18000,1);
    
    % Rescaling x and y;
    % inch to cm
    unit_c = 2.54;
    % pixel to cm
    alpha = unit_c*2400/(2048*100);                         %scaling number
    
    F1_X = A(1:18000,3);
    F1_Y = A(1:18000,4);
    
    % replace nan values
    F1_X = Replace_nan(F1_X);
    F1_Y = Replace_nan(F1_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F1_X = filtfilt(b,a,F1_X);
    F1_Y = filtfilt(b,a,F1_Y);
    
    F1_X = (alpha*F1_X) ;
    F1_Y = (alpha*F1_Y) ;
    
    minValue_F1_X = min(F1_X);
    minValue_F1_Y = min(F1_Y);
    maxValue_F1_X = max(F1_X);
    maxValue_F1_Y = max(F1_Y);
    beta1 = -((maxValue_F1_X+minValue_F1_X)/2);             %scaling number
    beta2 = -((maxValue_F1_Y+minValue_F1_Y)/2);    
    
    F1_X = (F1_X + beta1) ;
    F1_Y = (F1_Y + beta2) ;
    
    F2_X = A(1:18000,5);
    F2_Y = A(1:18000,6);
    
    % replace nan values
    F2_X = Replace_nan(F2_X);
    F2_Y = Replace_nan(F2_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F2_X = filtfilt(b,a,F2_X);
    F2_Y = filtfilt(b,a,F2_Y);
    
    F2_X = (alpha*F2_X) ;
    F2_Y = (alpha*F2_Y) ;
    
    minValue_F2_X = min(F2_X);
    minValue_F2_Y = min(F2_Y);
    maxValue_F2_X = max(F2_X);
    maxValue_F2_Y = max(F2_Y);
    beta3 = -((maxValue_F2_X+minValue_F2_X)/2);             %scaling number
    beta4 = -((maxValue_F2_Y+minValue_F2_Y)/2); 
    
    F2_X = (F2_X + beta3) ;
    F2_Y = (F2_Y + beta4) ;
    
    PF_D{i}.X1 = F1_X;
    PF_D{i}.Y1 = F1_Y;
    PF_D{i}.X2 = F2_X;
    PF_D{i}.Y2 = F2_Y;
end
save('Pairs_Fish_Dark.mat',"PF_D")

% Terns
for i = 1:n1
    filename = sprintf('Three_L_%d.csv', i);
    A = readmatrix(filename);
    time = A(1:18000,1);
    
    % Rescaling x and y;
    % inch to cm
    unit_c = 2.54;
    % pixel to cm
    alpha = unit_c*2400/(2048*100);                         %scaling number
    
    F1_X = A(1:18000,3);
    F1_Y = A(1:18000,4);
    
    % replace nan values
    F1_X = Replace_nan(F1_X);
    F1_Y = Replace_nan(F1_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F1_X = filtfilt(b,a,F1_X);
    F1_Y = filtfilt(b,a,F1_Y);
    
    F1_X = (alpha*F1_X) ;
    F1_Y = (alpha*F1_Y) ;
    
    minValue_F1_X = min(F1_X);
    minValue_F1_Y = min(F1_Y);
    maxValue_F1_X = max(F1_X);
    maxValue_F1_Y = max(F1_Y);
    beta1 = -((maxValue_F1_X+minValue_F1_X)/2);             %scaling number
    beta2 = -((maxValue_F1_Y+minValue_F1_Y)/2);    
    
    F1_X = (F1_X + beta1) ;
    F1_Y = (F1_Y + beta2) ;
    
    F2_X = A(1:18000,5);
    F2_Y = A(1:18000,6);
    
    % replace nan values
    F2_X = Replace_nan(F2_X);
    F2_Y = Replace_nan(F2_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F2_X = filtfilt(b,a,F2_X);
    F2_Y = filtfilt(b,a,F2_Y);
    
    F2_X = (alpha*F2_X) ;
    F2_Y = (alpha*F2_Y) ;
    
    minValue_F2_X = min(F2_X);
    minValue_F2_Y = min(F2_Y);
    maxValue_F2_X = max(F2_X);
    maxValue_F2_Y = max(F2_Y);
    beta3 = -((maxValue_F2_X+minValue_F2_X)/2);             %scaling number
    beta4 = -((maxValue_F2_Y+minValue_F2_Y)/2); 
    
    F2_X = (F2_X + beta3) ;
    F2_Y = (F2_Y + beta4) ;
    
    F3_X = A(1:18000,7);
    F3_Y = A(1:18000,8);
    
    % replace nan values
    F3_X = Replace_nan(F3_X);
    F3_Y = Replace_nan(F3_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F3_X = filtfilt(b,a,F3_X);
    F3_Y = filtfilt(b,a,F3_Y);
    
    F3_X = (alpha*F3_X) ;
    F3_Y = (alpha*F3_Y) ;
    
    minValue_F3_X = min(F3_X);
    minValue_F3_Y = min(F3_Y);
    maxValue_F3_X = max(F3_X);
    maxValue_F3_Y = max(F3_Y);
    beta5 = -((maxValue_F3_X+minValue_F3_X)/2);             %scaling number
    beta6 = -((maxValue_F3_Y+minValue_F3_Y)/2); 
    
    F3_X = (F3_X + beta5) ;
    F3_Y = (F3_Y + beta6) ;
    
    TF_B{i}.X1 = F1_X;
    TF_B{i}.Y1 = F1_Y;
    TF_B{i}.X2 = F2_X;
    TF_B{i}.Y2 = F2_Y;
    TF_B{i}.X3 = F3_X;
    TF_B{i}.Y3 = F3_Y;
end
save('Terns_Fish_Bright.mat',"TF_B")

for i = 1:n2
    filename = sprintf('Three_D_%d.csv', i);
    A = readmatrix(filename);
    time = A(1:18000,1);
    
    % Rescaling x and y;
    % inch to cm
    unit_c = 2.54;
    % pixel to cm
    alpha = unit_c*2400/(2048*100);                         %scaling number
    
    F1_X = A(1:18000,3);
    F1_Y = A(1:18000,4);
    
    % replace nan values
    F1_X = Replace_nan(F1_X);
    F1_Y = Replace_nan(F1_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F1_X = filtfilt(b,a,F1_X);
    F1_Y = filtfilt(b,a,F1_Y);
    
    F1_X = (alpha*F1_X) ;
    F1_Y = (alpha*F1_Y) ;
    
    minValue_F1_X = min(F1_X);
    minValue_F1_Y = min(F1_Y);
    maxValue_F1_X = max(F1_X);
    maxValue_F1_Y = max(F1_Y);
    beta1 = -((maxValue_F1_X+minValue_F1_X)/2);             %scaling number
    beta2 = -((maxValue_F1_Y+minValue_F1_Y)/2);    
    
    F1_X = (F1_X + beta1) ;
    F1_Y = (F1_Y + beta2) ;
    
    F2_X = A(1:18000,5);
    F2_Y = A(1:18000,6);
    
    % replace nan values
    F2_X = Replace_nan(F2_X);
    F2_Y = Replace_nan(F2_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F2_X = filtfilt(b,a,F2_X);
    F2_Y = filtfilt(b,a,F2_Y);
    
    F2_X = (alpha*F2_X) ;
    F2_Y = (alpha*F2_Y) ;
    
    minValue_F2_X = min(F2_X);
    minValue_F2_Y = min(F2_Y);
    maxValue_F2_X = max(F2_X);
    maxValue_F2_Y = max(F2_Y);
    beta3 = -((maxValue_F2_X+minValue_F2_X)/2);             %scaling number
    beta4 = -((maxValue_F2_Y+minValue_F2_Y)/2); 
    
    F2_X = (F2_X + beta3) ;
    F2_Y = (F2_Y + beta4) ;
    
    F3_X = A(1:18000,7);
    F3_Y = A(1:18000,8);
    
    % replace nan values
    F3_X = Replace_nan(F3_X);
    F3_Y = Replace_nan(F3_Y);
    
    % filter
    order = 4;
    wc = 0.4;
    [b,a] = butter(order, wc, 'low');
    F3_X = filtfilt(b,a,F3_X);
    F3_Y = filtfilt(b,a,F3_Y);
    
    F3_X = (alpha*F3_X) ;
    F3_Y = (alpha*F3_Y) ;
    
    minValue_F3_X = min(F3_X);
    minValue_F3_Y = min(F3_Y);
    maxValue_F3_X = max(F3_X);
    maxValue_F3_Y = max(F3_Y);
    beta5 = -((maxValue_F3_X+minValue_F3_X)/2);             %scaling number
    beta6 = -((maxValue_F3_Y+minValue_F3_Y)/2); 
    
    F3_X = (F3_X + beta5) ;
    F3_Y = (F3_Y + beta6) ;
    
    TF_D{i}.X1 = F1_X;
    TF_D{i}.Y1 = F1_Y;
    TF_D{i}.X2 = F2_X;
    TF_D{i}.Y2 = F2_Y;
    TF_D{i}.X3 = F3_X;
    TF_D{i}.Y3 = F3_Y;
end
save('Terns_Fish_Dark.mat',"TF_D")


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