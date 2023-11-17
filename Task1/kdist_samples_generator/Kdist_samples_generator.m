%% Uncorrelated K-distribution Samples generator
%Kdist_samples_generator.m

%This program generates Kdistributed samples with a given shape parameter

clear all
close all
clc


Data_length = 1e5;

Gaussian_data = randn(1,Data_length);

v = 0.3;
mu = 1;

Q = 0.5*(erfc(Gaussian_data/(2^0.5)));

Gamma_data = gaminv(Q,v,mu);

modulation = Gamma_data.* mu;

r = rand(Data_length,1);
kdata = -modulation .* log(r');

Kdist_data = kdata .* mu;

figure
%plot(Kdist_data(1:1000)')
histogram(Kdist_data)
title(['Kdistributed Samples with Shape parameter = ', num2str(v)])

%Estimated Samples Shape parameter using the samples Moments
Kdist_moments = Wa_NormalisedIntensityMoments(Kdist_data);
Est_Kdist_Shape_mo = (Kdist_moments(2)/2-1)^-1

%%
% Test fit simulated data with Kdist fitting methods - Using Sum Square
% Difference Fit

Kdist_data_norm = Kdist_data/mean(Kdist_data);

Kdist_data_norm_log = 10*log10(Kdist_data_norm);

Kdist_data_sort = sort(Kdist_data_norm_log);
Pfa_X= sort(Kdist_data_norm_log);

CDF_y = (1:length(Kdist_data_sort))./length(Kdist_data_sort);

Pfa_y = 1-CDF_y;

Pfa_Y = log10(Pfa_y);

pfa_pos_index = find(Pfa_X > 0);

Pfa_X_pos = Pfa_X(pfa_pos_index);
Pfa_Y_pos = Pfa_Y(pfa_pos_index);

[Pfa_numerical_array,nu_array] = RDV_kdist6(Pfa_X_pos);

RDV_overlay_plots=1;
RDV_pfa_seperate_Cmp=0;
[KdistShape,total_sum_sqr_diff] = RDV_fit_05(Pfa_X_pos,Pfa_numerical_array,Pfa_Y_pos,nu_array,RDV_overlay_plots,RDV_pfa_seperate_Cmp);

Kdist_moments = Wa_NormalisedIntensityMoments(Kdist_data);
Est_Kdist_Shape_mo = (Kdist_moments(2)/2-1)^-1


disp(['Desired Shape parameter of Kdist distribution ', num2str(v)])
disp(['Estimated Shape parameter from MSSD fitting ', num2str(KdistShape)])
