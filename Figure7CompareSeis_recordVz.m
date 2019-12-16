% first run Figure11aStaggeredFDTra.m,figure11bStaggeredFD.m,figure11cKspace to get seismic records,
% then run figure11compareSeimogramsVz.m, Figure11CompareSeis_recordVz.m; figure12compareSeimogramsTxx.m,figure12CompareSeis_recordTxx
% to get the figures in figure 11 and figure 12.
% This is only for the convenient of the reviewers.
clear;

clc
close all


load('StaggeredGridFDTra3.mat')
% figure;imagesc(-seis_recordVz(1:4000,45:end-45),[-0.0000002 0.0000002]);
% colormap gray
% xlabel('x/dx')
% ylabel('Travel time(ms)')
% title('')
plotimage(-seis_recordVz(1:4000,45:end-45))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')



load('StaggeredGridFD4.mat')
% figure;imagesc(-seis_recordVz(1:4000,45:end-45),[-0.0002 0.0002]);
% colormap gray
% xlabel('x/dx')
% ylabel('Travel time(ms)')
% title('')
plotimage(-seis_recordVz(1:4000,45:end-45))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

load('KSpaceLayer2.mat')
% figure;imagesc(-sign(real(seis_recordVz(1:4000,45:end-45))).*sqrt(real(seis_recordVz(1:4000,45:end-45)).^2+ imag(seis_recordVz(1:4000,45:end-45)).^2 ), [-0.0000002 0.0000002]);
% colormap gray
% xlabel('x/dx')
% ylabel('Travel time(ms)')

plotimage(-sign(real(seis_recordVz(1:4000,45:end-45))).*sqrt(real(seis_recordVz(1:4000,45:end-45)).^2+ imag(seis_recordVz(1:4000,45:end-45)).^2 ))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')