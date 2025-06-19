%% Computing magnitude and phase vs frequency in terms of fractions of pi

% Parameters
fs = 2000;              % Sampling frequency (Hz)
N = 4096;               % Number of samples
t = (0:N-1)/fs;         % Time vector

% Signal: two sine waves + noise
x = 16*sin(2*pi*1*t); % + 1.5*sin(2*pi*250*t) + 0.3*randn(1, N);

% Apply Hann window
% w = hann(N)'; % Uncomment when there is a discrete smpling like turbulent flow
w = ones(1,N);  % Uncomment if there is continuous sampling like sin wave
x_win = x .* w;

% Correct amplitude for energy loss due to window
A_corr = 1 / (sum(w)/N);  % amplitude correction factor

% FFT
X = fft(x_win);  % Compute FFT of the windowed signal
X_mag = abs(X) / N * A_corr;  % Normalize the magnitude and apply correction

% One-sided amplitude spectrum
X_amp = X_mag(1:N/2+1);
X_amp(2:end-1) = 2 * X_amp(2:end-1);  % Double the non-DC and non-Nyquist components

X_phase = zeros(1, N/2+1);  % Initialize phase array
for k = 1:N/2+1
    Re_Xk = real(X(k));  % Real part of X(k)
    Im_Xk = imag(X(k));  % Imaginary part of X(k)
    
    % Calculate phase using atan2 for correct quadrant
    X_phase(k) = atan2(Im_Xk, Re_Xk);  % Phase calculation
end

% Reconstruct the signal using amplitude and phase components
x_reconstructed = zeros(1, N);
for k = 1:N/2+1
    % Frequency for each bin
    f_k = (k-1) * fs / N;  % Frequency for the k-th bin
    % Reconstruct the signal for this frequency component using cosine
    x_reconstructed = x_reconstructed + X_amp(k) * cos(2 * pi * f_k * t + X_phase(k));  % Use cosine for reconstruction
end

% Frequency axis
f = (0:N/2) * fs / N;

% Plot the original signal and the reconstructed signal
figure;

% Plot original signal (before windowing)
subplot(3,1,1);
plot(t, x);  % Plot the original signal
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal');
grid on;

% Plot reconstructed signal (using amplitude and phase)
subplot(3,1,2);
plot(t, x_reconstructed);  % Plot the reconstructed signal
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal using Amplitude and Phase');
grid on;

subplot(3,1,3);
plot(t, X_amp(3) * cos(2 * pi * (3-1) * fs / N * t + X_phase(3)));  % Plot the reconstructed signal
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal using one frequency');
grid on;

% Plot magnitude spectrum
figure;
subplot(2,1,1);
plot(f, X_amp);  % Plot the magnitude spectrum
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('One-Sided Amplitude Spectrum with Hann Window');
grid on;

% Plot phase spectrum in terms of fractions of pi
subplot(2,1,2);
plot(f, X_phase);  % Plot the phase spectrum
xlabel('Frequency (Hz)');
ylabel('Phase');
title('Phase Spectrum (Unwrapped) in Fractions of \pi');
grid on;



%% Using Manual approach power spectrum density
% Set random seed for reproducibility
rng(42);  % You can use any integer

% Parameters
fs = 1000;              % Sampling frequency (Hz)
N = 2048;               % Total number of samples
t = (0:N-1)/fs;         % Time vector

% Signal: non-bin-aligned sinusoids + white noise
x = sin(2*pi*97.3*t) + 0.5*sin(2*pi*223.7*t) + 0.3*randn(1, N);  % Seeded noise

% Windowing setup
win_len = 512;                  % Window length
hop_size = 256;                 % 50% overlap
hann_win = hann(win_len)';      % Hann window (row vector)
U = sum(hann_win.^2) / win_len; % Window power normalization factor

% Segment the signal and compute FFT per segment
num_segments = floor((N - win_len)/hop_size) + 1;
nfft = win_len;
psd_accum = zeros(1, nfft);

for i = 1:num_segments
    idx_start = (i-1)*hop_size + 1;
    idx_end = idx_start + win_len - 1;
    x_seg = x(idx_start:idx_end) .* hann_win;
    
    X = fft(x_seg, nfft);
%     psd_seg = (abs(X).^2) / (fs * win_len * U);
    psd_seg = (X.*conj(X)) / (fs * win_len * U);
    psd_accum = psd_accum + psd_seg;
end

% Average PSD over segments
psd_avg = psd_accum / num_segments;

% One-sided PSD
psd_one_sided = psd_avg(1:nfft/2+1);
psd_one_sided(2:end-1) = 2 * psd_one_sided(2:end-1);

% Frequency axis
f = (0:nfft/2) * fs / nfft;

% Plot
figure;
plot(f, log(psd_one_sided));
xlabel('Frequency (Hz)');
ylabel('PSD (power/Hz)');
title('Manual One-Sided PSD with Hann Window (Seeded Noise)');
grid on;

%% Using pwelch
rng(42);
% Parameters
fs = 1000;           % Sampling frequency (Hz)
N = 2048;            % Number of samples
t = (0:N-1)/fs;      % Time vector

% Create signal: combination of sinusoids + wideband noise
x = sin(2*pi*97.3*t) + 0.5*sin(2*pi*223.7*t) + 0.3*randn(1, N);

% Use Welch method with Hann window to get smooth PSD
window = hann(512);
% window = ones(512,1); % You can try this and find the effect of energy lekeage
noverlap = 256;
nfft = 512;

[pxx, f] = pwelch(x, window, noverlap, nfft, fs);

% Plot
figure;
plot(f, log(pxx));  % PSD in dB/Hz
xlabel('Frequency (Hz)');
ylabel('PSD (power/Hz)');
title('One-Sided PSD of Signal with Continuous Spectrum');
grid on;

%% Using Manual approach for cross power spectrum density

% Set random seed for reproducibility
rng(123);

% Parameters
fs = 1000;              % Sampling frequency (Hz)
N = 2048;               % Number of samples
t = (0:N-1)/fs;         % Time vector

% Generate two correlated random signals (shared freq + noise)
x = sin(2*pi*120.7*t) + 0.5*randn(1, N);
y = 0.8*sin(2*pi*120.7*t + pi/4) + 0.5*randn(1, N);  % Same freq, phase shifted
% rho = 0.0;
% y = rho*x+sqrt(1-rho^2)*randn(1,N);

% Windowing setup
win_len = 512;
hop_size = 256;
hann_win = hann(win_len)';
U = sum(hann_win.^2) / win_len;

% Segment the signals and compute cross-spectrum
num_segments = floor((N - win_len)/hop_size) + 1;
nfft = win_len;
Sxy_accum = zeros(1, nfft);

for i = 1:num_segments
    idx_start = (i-1)*hop_size + 1;
    idx_end = idx_start + win_len - 1;
    
    x_seg = x(idx_start:idx_end) .* hann_win;
    y_seg = y(idx_start:idx_end) .* hann_win;
    
    X = fft(x_seg, nfft);
    Y = fft(y_seg, nfft);
    
    Sxy = (X .* conj(Y)) / (fs * win_len * U);  % Cross spectrum
    Sxy_accum = Sxy_accum + Sxy;
end

% Average over segments
Sxy_avg = Sxy_accum / num_segments;

% One-sided Cross Power Spectral Density
Sxy_one_sided = Sxy_avg(1:nfft/2+1);
Sxy_one_sided(2:end-1) = 2 * Sxy_one_sided(2:end-1);

% Frequency axis
f = (0:nfft/2) * fs / nfft;

% Plot magnitude of CPSD
figure;
plot(f, abs(Sxy_one_sided));
xlabel('Frequency (Hz)');
ylabel('|CPSD|');
title('Cross Power Spectral Density (Magnitude)');
grid on;

% Optional: Plot phase of CPSD
figure;
plot(f, angle(Sxy_one_sided));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase of Cross Power Spectral Density');
grid on;


%% Using cspd

% Set random seed for reproducibility
rng(123);

% Parameters
fs = 1000;              % Sampling frequency (Hz)
N = 2048;               % Number of samples
t = (0:N-1)/fs;         % Time vector

% Generate two correlated signals (shared sinusoid + noise)
x = sin(2*pi*120.7*t) ;%+ 0.5*randn(1, N);
y = 0.8*sin(2*pi*120.7*t + pi/4) ;%+ 0.5*randn(1, N);  % same freq, phase shifted
% rho = 0.1;
% y = rho*x+sqrt(1-rho^2)*randn(1,N);

% CPSD Parameters
window = hann(512);
nfft = 1024;
overlap = 256;

% Compute CPSD
[Sxy, f] = cpsd(x, y, window, overlap, nfft, fs);

% Plot magnitude of CPSD
figure;
plot(f, abs(Sxy));
xlabel('Frequency (Hz)');
ylabel('|CPSD|');
title('Cross Power Spectral Density (Magnitude) using cpsd');
grid on;

% Plot phase of CPSD
figure;
plot(f, angle(Sxy));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase of Cross Power Spectral Density using cpsd');
grid on;
