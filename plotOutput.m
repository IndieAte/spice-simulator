clc; clear all; close all;

T = readtable('output.csv');

frequency = T{:, 1};
amplitude = T{:, 2};
phase = T{:, 3};

figure;

title('Transfer Function')

yyaxis left;
semilogx(frequency, amplitude);
xlabel('Frequency / Hz');
ylabel('Amplitude / dB');

yyaxis right;
semilogx(frequency, phase, 'LineStyle', '--');
ylabel('Phase / Degrees');
grid on;
grid minor;