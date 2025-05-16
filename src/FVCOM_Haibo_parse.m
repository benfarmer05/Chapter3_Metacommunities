%% Script to parse Haibo's FVCOM data
%   Data was received 18 June 2024 through transfer to LSU HPC cluster, on
%   /project/holstein/haiboxu
% 18 June 2024

% Useful FVCOM resources:
%   - https://github.com/BigelowLab/fvcom
%   - https://github.com/moflaher/pyticle_tracker

clear;clc

cd /Volumes/UVI_Hydro_2019-2020/FVCOM_2019_Haibo-Xu/haiboxu

FVCOM = 'PRVI_20190102_0001.nc';

ncdisp(FVCOM)

partition = ncread(FVCOM, 'partition');