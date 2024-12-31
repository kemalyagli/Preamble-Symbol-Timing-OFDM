clear all; close all; clc 

%load('awgn_rayleigh_64_16.mat')
%load('awgn_rayleigh_128_16.mat')
load('awgn_rayleigh_256_32.mat')

figure;semilogy(EbN0dB, schmidl_mean_square_awgn,'o-.','LineWidth',2);hold on
semilogy(EbN0dB, minn_mean_square_awgn,'+-.','LineWidth',2);hold on
semilogy(EbN0dB, park_mean_square_awgn,'*-.','LineWidth',2);hold on
semilogy(EbN0dB, kim_mean_square_awgn,'x-.','LineWidth',2);hold on
semilogy(EbN0dB, ren_mean_square_awgn,'square-.','LineWidth',2);hold on
semilogy(EbN0dB, deneme_mean_square_awgn,'k-.','LineWidth',2);hold on
xlabel('SNR (dB)'); 
ylabel('Mean Square Error (MSE)'); 
legend("Schmidl's Method","Minn's Method", "Park's Method","Kim's Method","Ren's Method", "Proposed Method")
% axis([0 20 10^2 10^4])
grid on
grid minor

figure;semilogy(EbN0dB, schmidl_mean_square_rayleigh,'o-.','LineWidth',2);hold on
semilogy(EbN0dB, minn_mean_square_rayleigh,'+-.','LineWidth',2);hold on
semilogy(EbN0dB, park_mean_square_rayleigh,'*-.','LineWidth',2);hold on
semilogy(EbN0dB, kim_mean_square_rayleigh,'x-.','LineWidth',2);hold on
semilogy(EbN0dB, ren_mean_square_rayleigh,'square-.','LineWidth',2);hold on
semilogy(EbN0dB, deneme_mean_square_rayleigh,'k-.','LineWidth',2);hold on
xlabel('SNR (dB)'); 
ylabel('Mean Square Error (MSE)'); 
legend("Schmidl's Method","Minn's Method", "Park's Method","Kim's Method","Ren's Method", "Proposed Method")
% axis([0 20 10^2 10^4])
grid on
grid minor