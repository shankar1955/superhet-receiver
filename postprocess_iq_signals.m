%% postprocess_iq_signals.m
% Load I/Q outputs from workspace and extract features for neural network input.

% Sampling rate (Hz)
fs = evalin('base','fs');

% Band names
bands = {'Band1','Band2','Band3','Band4'};

% Preallocate feature table
features = struct();

for k = 1:numel(bands)
    b = bands{k};
    
    % Retrieve I and Q from workspace
  
    I = evalin('base', sprintf('out.%s_I',  b));
    Q = evalin('base', sprintf('out.%s_Q',  b));
    
    % Recombine complex baseband
    IQ = I + 1j*Q;
    
    % Time vector
    N = numel(IQ);
    t = (0:N-1)'/fs;
    
    % Inside your loop, replace the PSD section with:

% Determine window length (use at most N)
winLen = min(256, N);            % you can choose 256 or any power-of-two ≤ N
noverlap = floor(winLen/2);      % 50% overlap
nfft = max(512, winLen);         % at least 512 points FFT

% Compute PSD
[Pxx,f] = pwelch(IQ, hamming(winLen), noverlap, nfft, fs, 'centered');
    
    % 2. Constellation points (downsample for speed)
    idx = 1:round(N/1000):N;
    const_pts = IQ(idx);
    
    % 3. Statistical features
    mean_I   = mean(real(IQ));
    mean_Q   = mean(imag(IQ));
    var_I    = var(real(IQ));
    var_Q    = var(imag(IQ));
    skew_I   = skewness(real(IQ));
    skew_Q   = skewness(imag(IQ));
    kurt_I   = kurtosis(real(IQ));
    kurt_Q   = kurtosis(imag(IQ));
    
    % 4. Spectral centroid (Hz)
    centroid = sum(f.*Pxx)/sum(Pxx);
    
    % Store in struct
    features.(b).IQ        = IQ;
    features.(b).PSD.freq  = f;
    features.(b).PSD.pow   = Pxx;
    features.(b).ConstPts  = const_pts;
    features.(b).Stats     = [mean_I, mean_Q, var_I, var_Q, skew_I, skew_Q, kurt_I, kurt_Q];
    features.(b).Centroid  = centroid;
end

% Optionally convert features to matrix for neural network:
% Each row corresponds to one band: [mean_I mean_Q var_I var_Q skew_I skew_Q kurt_I kurt_Q centroid]
featureMatrix = zeros(numel(bands), 9);
for k = 1:numel(bands)
    b = bands{k};
    featureMatrix(k,1:8) = features.(b).Stats;
    featureMatrix(k,9)   = features.(b).Centroid;
end

% Plot I/Q Signals for all 4 bands using out.BandN_I and out.BandN_Q
bands = {'Band1', 'Band2', 'Band3', 'Band4'};
figure;

for k = 1:length(bands)
    b = bands{k};
    % Retrieve I and Q from workspace
    I = evalin('base', sprintf('out.%s_I', b));
    Q = evalin('base', sprintf('out.%s_Q', b));
    
    % Plot
    subplot(4, 1, k);
    plot(real(I), 'b'); hold on;
    plot(real(Q), 'r');
    grid on;
    title([b ' I (blue) and Q (red) Signals']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    legend('I', 'Q');
end


% Save features
assignin('base','features', features);
assignin('base','featureMatrix', featureMatrix);

fprintf('Post-processing complete. Features ready for neural network input.\n');

