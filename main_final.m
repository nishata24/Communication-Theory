% Seyun Kim, Lucia Rhode, Nishat Ahmed
% ECE300 Communication Theory Project 1: Error control coding

% 1. Rate 1/2 convolutional code
clc;
clear;
close all;

fprintf("Running 1/2 convolutional code\n")

% Define rate 1/2 convolutional encoder 
    % g0 = [1 0 1];
    % g1 = [1 1 1];

% Define parameters
%input = [0 0 1 1 0];
out = 2;
in = 1;
input = randi([0 1], 1, in*10^5);
contraintLength = 3;
probabilities = 0:0.01:0.5;
num = length(probabilities);
% encoded_length = length(input) * out
% tdepth formula taken from matlab
tdepth = 5*(contraintLength - 1);

% Initialize values
c_bsc = zeros(num, length(input) * out);
decoded = zeros(num, length(input));
bit_error_rate_1_2 = zeros(1, num);

% Define trellis
trellis = poly2trellis(contraintLength, [5 7]);

% Pass input to encoder
c = convenc(input, trellis);

for i = 1:length(probabilities)
    % Define binary symmetric channel with varying probabilities
    % Pass the encoded to each BSE
    c_bsc(i,:) =  bsc(c, probabilities(i));
 
    % Decode the output of BSE using Viterbi Algorithm
    decoded(i,:) = vitdec(c_bsc(i,:), trellis, tdepth, 'trunc', 'hard');

    % Compute bit error rate for each probability for BSE
    [number, ratio] = biterr(input, decoded(i,:));
    bit_error_rate_1_2(i) = ratio;
end

% Plot probability of BSC vs. BER
figure()
plot(probabilities, bit_error_rate_1_2, 'b')
title("Plot of probability of BSC vs BER of rate 1/2 Conv Code")
xlabel("Error Probability")
ylabel("Bit Error Rate (BER)")

% 2. Rate 2/3 Convolutional Code

fprintf("Running 2/3 convolutional code\n")

% Define rate 2/3 convolutional encoder 
    % g0 = [1 0 0 1 1; 0 0 0 0 0];
    % g1 = [1 1 1 0 1; 0 1 1 0 1];
    % g2 = [0 0 0 0 0; 0 1 0 1 1];

% Define parameters
% input = [0 0 1 1 1 0];
out = 3;
in = 2;
input = randi([0 1], 1, in*10^5);

contraintLength = [5 4];
% tdepth formula taken from matlab
% tdepth = round(7.5*(contraintLength  -1));

probabilities = 0:0.01:0.5;
num = length(probabilities);

% Initialize values
% c_bsc = zeros(num, length(input) / 2 * out);
% decoded = zeros(num, length(input));
bit_error_rate_2_3 = zeros(1, num);

% Define trellis
trellis = poly2trellis(contraintLength, [23 35 0; 0 5 13]);

% Pass input to encoder
c = convenc(input, trellis);

for i = 1:length(probabilities)
    % Define binary symmetric channel with varying probabilities
    % Pass the encoded to each BSE
    c_bsc =  bsc(c, probabilities(i));
 
    % Decode the output of BSE using Viterbi Algorithm
    decoded = vitdec(c_bsc, trellis, tdepth, 'trunc', 'hard');

    % Compute bit error rate for each probability for BSE
    [number, ratio] = biterr(input, decoded);
    bit_error_rate_2_3(i) = ratio;
end

% Plot probability of BSC vs. BER
figure()
plot(probabilities, bit_error_rate_2_3, 'b')
title("Plot of probability of BSC vs BER of rate 2/3 Conv Code")
xlabel("Error Probability")
ylabel("Bit Error Rate (BER)")


% 7-4 bch

fprintf("Running 7-4 BCH\n")

n = 7; %bits per codeword
k = 4; %bits per message
error_rate_7_4 = zeros(1, 51); %preallocate error rate array
probabilities = linspace(0, .5, 51); %range of error probabilities
bits = randi([0 1], 10^5, k); %random message bits generator
bch_msgTx = gf(bits); %message to encode as a Galois field array
for i = 1:51 
    enc_bch_msg = bchenc(bch_msgTx, n, k); %encode message using (n,k) bch encoder
    noisy_bch = bsc(enc_bch_msg, probabilities(i)); %pass encoded message through bsc
    bch_msgRx = bchdec(noisy_bch, n, k); %decode noisy code received from bsc using (n,k) bch decoder
    [number, BER] = biterr(bch_msgTx.x, bch_msgRx.x); %calculate number of errors and BER
    error_rate_7_4(i) = BER; %store BER values
    
end

% Plot BER vs. Probability
figure()
plot(probabilities, error_rate_7_4)
xlabel('Error Probability')
ylabel('Bit Error Rate (BER)')
title('BER vs. Probability of rate 7-4 BCH code')

% 15-7 bch
fprintf("Running 15-7 BCH\n")

n = 15; %bits per codeword
k = 7; %bits per message
error_rate_15_7 = zeros(1, 51); %preallocate error rate array
probabilities = linspace(0, .5, 51); %range of error probabilities
bits = randi([0 1], 10^5, k); %random message bits generator
bch_msgTx = gf(bits); %message to encode as a Galois field array
for i = 1:51 
    enc_bch_msg = bchenc(bch_msgTx, n, k); %encode message using (n,k) bch encoder
    noisy_bch = bsc(enc_bch_msg, probabilities(i)); %pass encoded message through bsc
    bch_msgRx = bchdec(noisy_bch, n, k); %decode noisy code received from bsc using (n,k) bch decoder
    [number, BER] = biterr(bch_msgTx.x, bch_msgRx.x); %calculate number of errors and BER
    error_rate_15_7(i) = BER; %store BER values
    
end 

% Plot BER vs. Probability
figure()
plot(probabilities, error_rate_15_7)
xlabel('Error Probability')
ylabel('Bit Error Rate (BER)')
title('BER vs. Probability of rate 15-7 BCH code')

% Plot of all four codes
figure()
hold on
plot(probabilities, bit_error_rate_1_2)
plot(probabilities, bit_error_rate_2_3)
plot(probabilities, error_rate_7_4)
plot(probabilities, error_rate_15_7)
xlabel('Error Probability')
ylabel('Bit Error Rate (BER)')
title('BER vs Probability of All Four Codes')
legend('rate 1/2 convolutional code', 'rate 2/3 convolutional code', 'rate 7-4 BCH code', 'rate 15-7 BCH code', 'Location','southeast')
