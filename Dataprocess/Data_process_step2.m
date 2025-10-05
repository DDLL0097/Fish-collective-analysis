clc
clear all
close all

dt = 1/30;

SF_D = load('Single_Fish_Dark.mat').SF_D;
SF_B = load('Single_Fish_Bright.mat').SF_B;

PF_D = load('Pairs_Fish_Dark.mat').PF_D;
PF_B = load('Pairs_Fish_Bright.mat').PF_B;

TF_D = load('Terns_Fish_Dark.mat').TF_D;
TF_B = load('Terns_Fish_Bright.mat').TF_B;

% Filtering
order = 4;
wc = 0.4;
[b, a] = butter(order, wc, 'low');

SF_B_X = cell(size(SF_B)); SF_B_Y = cell(size(SF_B));
SF_D_X = cell(size(SF_D)); SF_D_Y = cell(size(SF_D));

PF_B_X1 = cell(size(PF_B)); PF_B_Y1 = cell(size(PF_B));
PF_B_X2 = cell(size(PF_B)); PF_B_Y2 = cell(size(PF_B));
PF_D_X1 = cell(size(PF_D)); PF_D_Y1 = cell(size(PF_D));
PF_D_X2 = cell(size(PF_D)); PF_D_Y2 = cell(size(PF_D));

TF_B_X1 = cell(size(TF_B)); TF_B_Y1 = cell(size(TF_B));
TF_B_X2 = cell(size(TF_B)); TF_B_Y2 = cell(size(TF_B));
TF_B_X3 = cell(size(TF_B)); TF_B_Y3 = cell(size(TF_B));
TF_D_X1 = cell(size(TF_D)); TF_D_Y1 = cell(size(TF_D));
TF_D_X2 = cell(size(TF_D)); TF_D_Y2 = cell(size(TF_D));
TF_D_X3 = cell(size(TF_D)); TF_D_Y3 = cell(size(TF_D));


%% =====================================================================
%% single fish analysis
% Bright
for i=1:length(SF_B)
    X = SF_B{i}.X;
    Y = SF_B{i}.Y;
    X_1L = X;
    Y_1L = Y;

    SF_B_X{i} = X_1L;
    SF_B_Y{i} = Y_1L;
end

% Dark
for i=1:length(SF_D)
    X = SF_D{i}.X;
    Y = SF_D{i}.Y;
    X_1D = X;
    Y_1D = Y;

    SF_D_X{i} = X_1D;
    SF_D_Y{i} = Y_1D;
end


%% =====================================================================
%% Pairs fish analysis
% Bright
for i=1:length(PF_B)
    
    X1 = PF_B{i}.X1; Y1 = PF_B{i}.Y1;
    X2 = PF_B{i}.X2; Y2 = PF_B{i}.Y2;

    X1_2L = X1; Y1_2L = Y1;
    X2_2L = X2; Y2_2L = Y2;

    PF_B_X1{i} = X1_2L;
    PF_B_Y1{i} = Y1_2L;
    PF_B_X2{i} = X2_2L;
    PF_B_Y2{i} = Y2_2L;
end

% Dark
for i=1:length(PF_D)
    
    X1 = PF_D{i}.X1; Y1 = PF_D{i}.Y1;
    X2 = PF_D{i}.X2; Y2 = PF_D{i}.Y2;

    X1_2D = X1; Y1_2D = Y1;
    X2_2D = X2; Y2_2D = Y2;

    PF_D_X1{i} = X1_2D;
    PF_D_Y1{i} = Y1_2D;
    PF_D_X2{i} = X2_2D;
    PF_D_Y2{i} = Y2_2D;
end


%% =====================================================================
%% Terns fish analysis
% Bright
for i=1:length(TF_B)
    
    X1 = TF_B{i}.X1; Y1 = TF_B{i}.Y1;
    X2 = TF_B{i}.X2; Y2 = TF_B{i}.Y2;
    X3 = TF_B{i}.X3; Y3 = TF_B{i}.Y3;
   
    X1_3L = X1; Y1_3L = Y1;
    X2_3L = X2; Y2_3L = Y2;
    X3_3L = X3; Y3_3L = Y3;

    TF_B_X1{i} = X1_3L;
    TF_B_Y1{i} = Y1_3L;
    TF_B_X2{i} = X2_3L;
    TF_B_Y2{i} = Y2_3L;
    TF_B_X3{i} = X3_3L;
    TF_B_Y3{i} = Y3_3L;
end

% Dark
for i=1:length(TF_D)
    
    X1 = TF_D{i}.X1; Y1 = TF_D{i}.Y1;
    X2 = TF_D{i}.X2; Y2 = TF_D{i}.Y2;
    X3 = TF_D{i}.X3; Y3 = TF_D{i}.Y3;

    X1_3D = X1; Y1_3D = Y1;
    X2_3D = X2; Y2_3D = Y2;
    X3_3D = X3; Y3_3D = Y3;

    TF_D_X1{i} = X1_3D;
    TF_D_Y1{i} = Y1_3D;
    TF_D_X2{i} = X2_3D;
    TF_D_Y2{i} = Y2_3D;
    TF_D_X3{i} = X3_3D;
    TF_D_Y3{i} = Y3_3D;
end

% save
save('Single_Fish_Bright_XY.mat', 'SF_B_X', 'SF_B_Y');
save('Single_Fish_Dark_XY.mat',   'SF_D_X', 'SF_D_Y');

save('Pairs_Fish_Bright_XY.mat',   'PF_B_X1', 'PF_B_Y1', 'PF_B_X2', 'PF_B_Y2');
save('Pairs_Fish_Dark_XY.mat',     'PF_D_X1', 'PF_D_Y1', 'PF_D_X2', 'PF_D_Y2');

save('Terns_Fish_Bright_XY.mat',   'TF_B_X1', 'TF_B_Y1', 'TF_B_X2', 'TF_B_Y2', 'TF_B_X3', 'TF_B_Y3');
save('Terns_Fish_Dark_XY.mat',     'TF_D_X1', 'TF_D_Y1', 'TF_D_X2', 'TF_D_Y2', 'TF_D_X3', 'TF_D_Y3');