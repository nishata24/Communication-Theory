%% ECE300 Project 5
%% Lucia Rhode, Nishat Ahmed, Seyun Kim

clc;
clear;
close all;

%% MIMO 2x2 Link 
Nt = 2; %Number of transmit antennas
Nr = 2; %Number of receive antennas

%Define flat fading gains
H1 = (1/sqrt(2))*(randn(Nr,Nt)+1i*randn(Nr,Nt));
H2 = (1/sqrt(2))*(randn(Nr,Nt)+1i*randn(Nr,Nt));
H3 = (1/sqrt(2))*(randn(Nr,Nt)+1i*randn(Nr,Nt));

%Define transmitted signal
N = 10^6;
a1 = randi([0 1], 1, N);
a2 = randi([0 1], 1, N);

%BPSK modulation order
M = 2;

%perform BPSK modulation
x1 = pskmod(a1, M); 
x2 = pskmod(a2, M);
x = [x1; x2];

%pass through channel and add noise 
n1 = wgn(1, N, 10);
n2 = wgn(1, N, 10);
n = [n1; n2];

%Pre-coding 

%obtain U, S, V
[U, S, V] = svd(H1);

%output
x_ = V\x;
n_ = U'*n;
y_ = S*x_+n_;

%receiver shaping
y = U'\y_;

%demodulate 
y1 = pskdemod(y(1, :), M);
y2 = pskdemod(y(2, :), M);
y = [y1; y2];
BER1 = biterr([a1; a2], y);

%data rate
signal_power = sum(abs(H1).^2, 'all');
noise_var = .1;
SNR = signal_power/noise_var;
R1 = log2(det(eye(2) + SNR * H1 * H1'));
%R1 = 9.7857 + 0.0000i

%Zero-Forcing 
%output 
y_z = H2*x+n;

%receive
y_z = H2\y_z;

%demodulation
y_z1 = pskdemod(y_z(1,:), M);
y_z2 = pskdemod(y_z(2,:), M);
y_z = [y_z1; y_z2];
BER2 = biterr([a1; a2], y_z);

%data rate
signal_power = sum(abs(H2).^2, 'all');
noise_var = .1;
SNR = signal_power/noise_var;
R2 = log2(det(eye(2) + SNR * H2' * H2));
%R2 = 8.5941 + 0.0000i

%MMSE
a = [a1; a2];
y_mm = H3*x + n;
length = 5;
weights = zeros(length+1,1);
for i = 1:size(y_mm,1)
    for j = length+1:size(y_mm,2)
        xn = flip(y_mm(i,j-length:j))';
        yn = xn'*weights;
        en = a(i,j) - yn;
        weights = weights + .001*en*xn;
    end
    y(i,:) = filter(weights, 1, y(i,:));
end

%demodulation
y_z1 = pskdemod(y_mm(1, :), M);
y_z2 = pskdemod(y_mm(2, :), M);
y_mm = [y_z1; y_z2];
BER3 = biterr(a, y_mm);

%data rate
signal_power = sum(abs(H3).^2, 'all');
noise_var = .1;
SNR = signal_power/noise_var;
W = inv(H3' * H3 + noise_var / SNR * eye(2)) * H3'; %MMSE filter 
R3 = log2(det(eye(2) + SNR * H3 * W * H3'));
%R3 = 10.0245 - 2.9108i

%MMSE receiver provides the greatest data rate. This is because MMSE can
%provide good performance at all SNR values, even though it has high
%complexity. On the contrary, pre-coding provides good results at high SNR
%values but can suffer from sensitivity to channel estimation errors.
%Similarly, ZF provides good results at moderate to high SNR values but can
%have poor results for low SNR values. 

%% OFDM 
num_subcarriers = 64;
num_subcarriers_tx = 48; 
mu = 16;
r1 = 1/2;
r2 = 2/3;
r3 = 3/4;
M = 2;
N2 = 1000;
a = randi([0 1],1,N2);
bpsk = pskmod(a,M);
bandwidth = 20e6; 
Ts = 1/bandwidth;
Bn = bandwidth/num_subcarriers; %312.5kHz
Tm_max = mu*Ts; %.8us
Tn = (num_subcarriers+mu)*Ts; %4us
rn = log2(M/Tn); 
rate_min = num_subcarriers_tx/Tn*r1; %6Mbps
rate_max = num_subcarriers_tx*r3*6/Tn; %54Mbps

%16 QAM modulator
M = 16;
bits = randi([0 M-1], 1, N2);
x_qam = qammod(bits, M);

%serial to parallel
x_qam = x_qam.';

%IFFT, cyclic prefix, and parallel to serial
x = ifft(x_qam, N2);
x_ = [x(N2-mu+1:N2); x].';

%frequency selective channel
h = fliplr(linspace(0.1, 1, 17));
H = zeros(N2, N2+mu); 
for i = 1:size(x_, 2)-mu
    H(i, i:i+M) = h;
end

noise = wgn(N2, 1, 0);

%received
y = H*flipud(x_.') + noise;
y = flipud(y);

%FFT
Y = fft(y, N2);

%parallel to serial
Y = Y.';

%Zero-forcing
Y = ifft(fft(Y, N2)./fft(h, N2), N2);
bits_rx = qamdemod(Y, M);
BER4 = biterr(bits, bits_rx);

%MMSE
Y = ifft(fft(Y, N2)./(fft(h, N2)+fft(noise.', N2)), N2);
bits_rx = qamdemod(Y, M);
BER5 = biterr(bits, bits_rx);









