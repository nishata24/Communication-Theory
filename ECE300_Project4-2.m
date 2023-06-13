%% ECE300 Project 4
%% Nishat Ahmed, Seyun Kim, Lucia Rhode

clc;
clear;
close all;

%independent Rayleigh gains
N= 10^6;
h0=(1/sqrt(2))*(randn(N,1)+1i*randn(N,1));
h1=(1/sqrt(2))*(randn(N,1)+1i*randn(N,1));

%BPSK modulation order
M = 2;

%transmitter 
bits = randi([0 1], N, 1); %generate random bits 


bpsk = pskmod(bits, M); %perform BPSK modulation 

%% no diversity (1 Tx, 1 Rx)

%preallocate error rate array
error_rate0 = zeros(1, 51); 
y_decision0 = zeros(N, 1);

count = 0;

for snr = 0:50
    y0 = awgn(h0.*bpsk, snr);
   
    y_combined = (conj(h0).*y0); 
 
    for i = 1:numel(y_combined)
        if y_combined(i) > 0
            y_decision0(i, 1) = 0;
        end
        if y_combined(i) < 0 
            y_decision0(i, 1) = 1;
        end
    end

    %bit error rate 
    BER = biterr(y_decision0, bits);
    count = count + 1;
    error_rate0(count) = BER;
end

%% MRRC (1Tx, 2Rx)

%preallocate error rate array
error_rate = zeros(1, 51); 
y_decision = zeros(N, 1);

count = 0;

for snr = 0:50
    y0 = awgn(h0.*bpsk, snr);
    y1 = awgn(h1.*bpsk, snr);
    y_combined = (conj(h0).*y0)+(conj(h1).*y1); 
 
    for i = 1:numel(y_combined)
        if y_combined(i) > 0
            y_decision(i, 1) = 0;
        end
        if y_combined(i) < 0 
            y_decision(i, 1) = 1;
        end
    end

    %bit error rate 
    BER = biterr(y_decision, bits);
    count = count + 1;
    error_rate(count) = BER;
end

%% MRRC (1Tx, 4Rx)

h2=(1/sqrt(2))*(randn(N,1)+1i*randn(N,1));
h3=(1/sqrt(2))*(randn(N,1)+1i*randn(N,1));

%preallocate error rate array
error_rate2 = zeros(1, 51); 
y_decision2 = zeros(N, 1);

count = 0;

for snr = 0:50
    y0 = awgn(h0.*bpsk, snr);
    y1 = awgn(h1.*bpsk, snr);
    y2 = awgn(h2.*bpsk, snr);
    y3 = awgn(h3.*bpsk, snr);

    y_combined = (conj(h0).*y0)+(conj(h1).*y1)+(conj(h2).*y2)+(conj(h3).*y3); 
 
    for i = 1:numel(y_combined)
        if y_combined(i) > 0
            y_decision2(i, 1) = 0;
        end
        if y_combined(i) < 0 
            y_decision2(i, 1) = 1;
        end
    end

    %bit error rate 
    BER = biterr(y_decision2, bits);
    count = count + 1;
    error_rate2(count) = BER;
end

%% New Scheme (2 Tx, 1Rx)

bpsk1 = pskmod(bits, M)/sqrt(2); %perform BPSK modulation with same transmitted power as MRCC
bits2 = randi([0 1], N, 1); %generate random bits 
bpsk2 = pskmod(bits2, M)/sqrt(2); %perform BPSK modulation with same transmitted power as MRCC

%preallocate error rate array
error_rate3 = zeros(1, 51); 
y_decision3 = zeros(N, 1);
y_decision4 = zeros(N, 1);

count = 0;

for snr = 0:50
    y0 = awgn((h0.*bpsk1)+(h1.*bpsk2), snr);
    y1 = awgn(-(h0.*conj(bpsk2))+(h1.*conj(bpsk1)), snr);
   
    y_combined = (conj(h0).*y0)+(h1.*conj(y1)); 
    y_combined2 = (conj(h1).*y0)-(h0.*conj(y1));

    for i = 1:numel(y_combined)
        if y_combined(i) > 0
            y_decision3(i, 1) = 0;
        end
        if y_combined(i) < 0 
            y_decision3(i, 1) = 1;
        end
    end

    for i = 1:numel(y_combined2)
        if y_combined2(i) > 0
            y_decision4(i, 1) = 0;
        end
        if y_combined2(i) < 0 
            y_decision4(i, 1) = 1;
        end
    end

    %bit error rate 
    BER = biterr([y_decision3; y_decision4], [bits; bits2]);
    count = count + 1;
    error_rate3(count) = BER;
end

%% New Scheme (2 Tx, 2 Rx)

%preallocate error rate array
error_rate4 = zeros(1, 51); 
y_decision5 = zeros(N, 1);
y_decision6 = zeros(N, 1);

count = 0;

for snr = 0:50
    y0 = awgn((h0.*bpsk1)+(h1.*bpsk2), snr);
    y1 = awgn(-(h0.*conj(bpsk2))+(h1.*conj(bpsk1)), snr);
    y2 = awgn((h2.*bpsk1)+(h3.*bpsk2), snr);
    y3 = awgn(-(h2.*conj(bpsk2))+(h3.*conj(bpsk1)), snr);
   
    y_combined = (conj(h0).*y0)+(h1.*conj(y1))+(conj(h2).*y2)+(h3.*conj(y3)); 
    y_combined2 = (conj(h1).*y0)-(h0.*conj(y1))+(conj(h3).*y2)-(h2.*conj(y3));

    for i = 1:numel(y_combined)
        if y_combined(i) > 0
            y_decision5(i, 1) = 0;
        end
        if y_combined(i) < 0 
            y_decision5(i, 1) = 1;
        end
    end

    for i = 1:numel(y_combined2)
        if y_combined2(i) > 0
            y_decision6(i, 1) = 0;
        end
        if y_combined2(i) < 0 
            y_decision6(i, 1) = 1;
        end
    end

    %bit error rate 
    BER = biterr([y_decision5; y_decision6], [bits; bits2]);
    count = count + 1;
    error_rate4(count) = BER;
end

snr = 0:50;
figure()
semilogy(snr, error_rate0/N);
hold on 
semilogy(snr, error_rate/N);
hold on
semilogy(snr, error_rate2/N);
hold on 
semilogy(snr, error_rate3/N);
hold on 
semilogy(snr, error_rate4/N);
hold off
xlabel('Signal to Noise Ratio (SNR)')
ylabel('Bit Error Rate (BER)')
title('BER vs. SNR')
legend('no diversity (1 Tx, 1 Rx)', 'MRCC (1 Tx, 2 Rx)', 'MRCC (1 Tx, 4 Rx)', 'new scheme (2 Tx, 1 Rx)', 'new scheme (2Tx, 2 Rx)')
