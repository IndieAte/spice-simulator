clc; clear all; close all;

T = readtable('output.csv');

frequency = T{:, 1};
amplitude = T{:, 2};
phase = T{:, 3};

figure;
semilogx(frequency, amplitude);
title('Amplitude versus Frequency');
xlabel('Frequency / Hz');
ylabel('Amplitude / dB');
grid on;
grid minor;

figure;
semilogx(frequency, phase);
title('Phase versus Frequency');
xlabel('Frequency / Hz');
ylabel('Phase / Degrees');
grid on;
grid minor;