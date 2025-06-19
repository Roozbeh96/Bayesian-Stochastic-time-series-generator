% Step 1: Generate two 1D signals
n = 256;
x = linspace(0, 2*pi, n);

% Original field (target spectrum): smooth signal
signal_A = randn(1, n);%sin(5*x) + 0.5*sin(10*x);

% Generated field: random noise
signal_G = randn(1, n);

% Step 2: Compute FFTs
fft_A = fft(signal_A);
fft_G = fft(signal_G);

% Get amplitudes
amp_A = abs(fft_A);
amp_G = abs(fft_G);

% Step 3: Replace amplitude of generated field with that of original field
epsilon = 1e-8;
phase_G = fft_G ./ (amp_G + epsilon);
fft_G_filtered = phase_G .* amp_A;

% Step 4: Inverse FFT to get filtered signal
signal_G_filtered = ifft(fft_G_filtered);

% Step 5: Plotting
figure;

subplot(3,1,1);
plot(x, signal_A, 'b', 'LineWidth', 1.5);
title('Original Signal (Target)');
xlabel('x'); ylabel('Amplitude');

subplot(3,1,2);
plot(x, signal_G, 'r', 'LineWidth', 1.5);
title('Generated Signal (Before Filtering)');
xlabel('x'); ylabel('Amplitude');

subplot(3,1,3);
plot(x, signal_G_filtered, 'g', 'LineWidth', 1.5);
title('Filtered Generated Signal (Matched Spectrum)');
xlabel('x'); ylabel('Amplitude');

% Optional: Compare Power Spectra
figure;
hold on;
plot(log10(amp_A.^2), 'b', 'DisplayName', 'Target Spectrum');
plot(log10(amp_G.^2), 'r', 'DisplayName', 'Before Filtering');
plot(log10(abs(fft(signal_G_filtered)).^2), 'g', 'DisplayName', 'After Filtering');
legend;
title('Log Power Spectrum Comparison');
xlabel('Frequency Index');
ylabel('Log Power');
