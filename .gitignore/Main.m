%% Load data sets and array information.
clc;clear;
load('data.mat');
load('ntk.mat');% array info, including channel no, locations, sampling rate, etc.
%% Set parameters for extraction
param.min_points=10;
param.jump_space=80;
param.jump_frame = 5;
param.init_delay=5;
%% Axon extraction
for i =1:length(neurons)
     [neurons{i}.branch] = axon_velocity_auto(neurons{i},ntk,param,1);
end
