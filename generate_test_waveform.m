function generate_test_waveform()
% Parameters
fs = 40e6; % Sampling rate: 40 MSPS
T = 1e-3; % Duration: 1 ms
t = (0:1/fs:T-1/fs)'; % Time vector (column)

% Carrier frequencies for each band
freqs = [400e6, 2.45e9, 3.5e9, 5.9e9];

% Initialize test signal
input_signal = zeros(size(t));

% Add cosine carriers
for f = freqs
    input_signal = input_signal + 0.1*cos(2*pi*f*t);
end

% ——— BROADBAND EMI INJECTION START ———
fprintf('=== BROADBAND EMI INJECTION FOR POWER SYSTEMS ===\n');

% EMI Parameters for different interference types
emi_types = {
    struct('name', 'Power_Line_Harmonics', 'center_freq', 0, 'bandwidth', 2000, 'power', -40), % 2 kHz BW around DC
    struct('name', 'Corona_Discharge', 'center_freq', 0, 'bandwidth', 100000, 'power', -50), % 100 kHz broadband
    struct('name', 'Switching_Transients', 'center_freq', 0, 'bandwidth', 10000000, 'power', -35) % 10 MHz broadband
};

signal_before_emi = input_signal;

for emi_idx = 1:length(emi_types)
    emi = emi_types{emi_idx};
    fprintf('Adding %s...\n', emi.name);
    
    % Generate broadband noise
    noise_signal = randn(size(t));
    
    % Design bandpass filter for EMI bandwidth
    if emi.bandwidth < fs/2  % Nyquist check
        if emi.center_freq == 0  % Low-pass for baseband EMI
            % Low-pass filter for power line and corona EMI
            cutoff_freq = min(emi.bandwidth/2, fs/2 - 1000); % Leave some margin
            [b_emi, a_emi] = butter(4, cutoff_freq/(fs/2), 'low');
        else
            % Bandpass filter (for future use with specific center frequencies)
            f_low = max(100, emi.center_freq - emi.bandwidth/2);
            f_high = min(fs/2-1000, emi.center_freq + emi.bandwidth/2);
            [b_emi, a_emi] = butter(4, [f_low f_high]/(fs/2), 'bandpass');
        end
        
        % Apply filter to noise
        filtered_noise = filtfilt(b_emi, a_emi, noise_signal);
        
        % Scale to desired power level
        power_scale = 10^(emi.power/20);  % Convert dB to linear scale
        scaled_emi = power_scale * filtered_noise;
        
        % Add to input signal
        input_signal = input_signal + scaled_emi;
        
        fprintf('  Bandwidth: %.1f kHz, Power: %.1f dB, RMS: %.6f\n', ...
                emi.bandwidth/1000, emi.power, rms(scaled_emi));
    else
        fprintf('  Skipped - bandwidth exceeds Nyquist limit\n');
    end
end

% Add impulsive switching noise (time-domain bursts)
fprintf('Adding impulsive switching events...\n');
num_bursts = 5;  % Number of switching events
for burst_idx = 1:num_bursts
    % Random burst timing
    burst_time = rand() * T;
    burst_start = round(burst_time * fs);
    burst_duration = round(10e-6 * fs);  % 10 microsecond bursts
    
    if burst_start + burst_duration <= length(t)
        % Create exponentially decaying burst
        burst_envelope = exp(-linspace(0, 5, burst_duration)');
        burst_noise = 0.05 * randn(burst_duration, 1) .* burst_envelope;
        
        % Add burst to signal
        burst_end = burst_start + burst_duration - 1;
        input_signal(burst_start:burst_end) = input_signal(burst_start:burst_end) + burst_noise;
        
        fprintf('  Burst %d at %.2f ms (duration: %.1f μs)\n', ...
                burst_idx, burst_time*1000, burst_duration/fs*1e6);
    end
end

% Add thermal noise floor
thermal_noise_power = -60; % dB
thermal_scale = 10^(thermal_noise_power/20);
thermal_noise = thermal_scale * randn(size(t));
input_signal = input_signal + thermal_noise;

fprintf('Added thermal noise floor at %.1f dB\n', thermal_noise_power);

% ——— EMI INJECTION END ———

% Calculate power changes
power_before = mean(abs(signal_before_emi).^2);
power_after = mean(abs(input_signal).^2);
power_increase = power_after - power_before;

fprintf('\nBROADBAND EMI INJECTION SUMMARY:\n');
fprintf('  Signal power BEFORE: %.6e\n', power_before);
fprintf('  Signal power AFTER:  %.6e\n', power_after);
fprintf('  Power increase:      %.6e (%.2f%%)\n', power_increase, 100*power_increase/power_before);

% Verify EMI injection with spectrum analysis
fprintf('\nEMI Injection Verification:\n');
N_check = length(input_signal);
f_check = (0:N_check-1)*(fs/N_check);
X_before = abs(fft(signal_before_emi));
X_after = abs(fft(input_signal));

% Check power increase in low frequencies (0-10 MHz)
low_freq_mask = f_check <= 10e6;
power_increase_low = sum(X_after(low_freq_mask).^2) - sum(X_before(low_freq_mask).^2);

fprintf('  Low-frequency power increase: %.2e\n', power_increase_low);
fprintf('  EMI successfully added: %s\n', power_increase_low > 1e-6);

% Save to workspace
input_signal_full = [t, input_signal];
assignin('base','input_signal_full',input_signal_full);
assignin('base','input_signal',input_signal);
assignin('base','time_vector',t);
assignin('base','fs',fs);
assignin('base','signal_before_emi', signal_before_emi);

% Create verification plot
figure('Name', 'Broadband EMI Injection Verification');

% Spectrum comparison
subplot(2,1,1);
f_plot_MHz = f_check/1e6;  % Convert to MHz
X_before_db = 20*log10(X_before + 1e-12);
X_after_db = 20*log10(X_after + 1e-12);

% Plot low-frequency range (0-20 MHz)
freq_mask_plot = f_plot_MHz <= 20;
plot(f_plot_MHz(freq_mask_plot), X_before_db(freq_mask_plot), 'b-', 'DisplayName', 'Before EMI');
hold on;
plot(f_plot_MHz(freq_mask_plot), X_after_db(freq_mask_plot), 'r-', 'DisplayName', 'After EMI');

xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('EMI Injection Verification (0-20 MHz)');
legend;
grid on;

% Time domain comparison
subplot(2,1,2);
t_plot_us = t(1:min(4000, length(t))) * 1e6;  % First 4000 samples in microseconds
plot(t_plot_us, signal_before_emi(1:length(t_plot_us)), 'b-', 'DisplayName', 'Before EMI');
hold on;
plot(t_plot_us, input_signal(1:length(t_plot_us)), 'r-', 'DisplayName', 'After EMI');

xlabel('Time (μs)');
ylabel('Amplitude');
title('Time Domain - EMI Effects');
legend;
grid on;

fprintf('Test waveform generated with broadband EMI:\n');
fprintf('  Duration: %.2f ms\n', T*1e3);
fprintf('  Samples: %d\n', length(t));
fprintf('  Sampling: %.1f MSPS\n', fs/1e6);
fprintf('  EMI types: %d different interference sources\n', length(emi_types));

input_signal_full = [t, input_signal];

end