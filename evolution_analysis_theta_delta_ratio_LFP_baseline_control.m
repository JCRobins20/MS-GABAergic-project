
% Example time stamps from events with code 128
laser_ts = cell2mat(root_stim.event(:, 4));
on_ts = cell2mat(root_stim.event(laser_ts == 128, 2));

% Calculate start and end times for off periods
off_ts = on_ts + 30; % Start of "off" period
off_end = off_ts + 60; % End of "off" period after an additional 30 seconds

% Ensure the channel number provided is valid
if channel < 1 || channel > length(root_base.b_lfp) || isnan(channel)
    error('Invalid channel number.');
end

% Channel and data extraction for actual experiment
data = root_stim.b_lfp(channel).signal;
ts = root_stim.b_lfp(channel).ts;

% Channel and data extraction for control
data1 = root_base.b_lfp(channel).signal;  % Assuming root_base is defined and loaded similarly
ts1 = root_base.b_lfp(channel).ts;

% Example setup -- ensure these are defined as per your actual data setup
Fs = 500;  % Sampling frequency
windowsize = 30;  % 30 seconds
windowinc = 2;  % Move 2 seconds at a time

% Parameters for spectral analysis, adjusted for smaller data segments
params.tapers = [3, 5];  % Lower time-bandwidth product and fewer tapers
params.Fs = Fs;
params.fpass = [4 12];  % Theta band
params.pad = 1;  % Padding for spectral analysis

% Generate 1000 random stimulation times for control
exp_start = root_base.epoch(1);
exp_end = root_base.epoch(2) - 90; % Ensures full window analysis for last timestamp
rand_on_ts = exp_start + (exp_end - exp_start) .* rand(1000, 1);
rand_on_ts = sort(rand_on_ts); % Sort times

% Calculate maximum number of windows to analyze
max_windows = floor((90 - windowsize) / windowinc) + 1;

% Initialize matrix for storing power values
theta_powers = NaN(max_windows, length(on_ts));
rand_theta_powers = NaN(max_windows, length(rand_on_ts));

% Analyze actual data
for j = 1:length(on_ts)
    for i = 0:max_windows-1
        window_start = on_ts(j) + i * windowinc;
        window_end = window_start + windowsize;
        if window_end <= max(ts)
            indices = ts >= window_start & ts < window_end;
            temp_sig = data(indices);
            if ~isempty(temp_sig)
                [S, f] = mtspectrumc(detrend(temp_sig), params);
                theta_band_power = mean(S(f >= 4 & f <= 12));
                theta_powers(i+1, j) = theta_band_power;
            end
        end
    end
end

% Analyze control data
for j = 1:length(rand_on_ts)
    for i = 0:max_windows-1
        window_start = rand_on_ts(j) + i * windowinc;
        window_end = window_start + windowsize;
        if window_end <= max(ts1)
            indices = ts1 >= window_start & ts1 < window_end;
            temp_sig = data1(indices);
            if ~isempty(temp_sig) && length(temp_sig) > (2 * max(params.tapers(1, 1:2)))  % Ensure sufficient length for spectral analysis
                [S, f] = mtspectrumc(detrend(temp_sig), params);
                rand_theta_band_power = mean(S(f >= 4 & f <= 12));
                rand_theta_powers(i+1, j) = rand_theta_band_power;
            end
        end
    end
end

% Calculate mean and SEM of powers across all on periods
mean_powers = nanmean(theta_powers, 2);
sem_powers = nanstd(theta_powers, 0, 2) ./ sqrt(sum(~isnan(theta_powers), 2));
% Calculate mean and SEM of powers across all on periods
rand_mean_powers = nanmean(rand_theta_powers, 2);
rand_sem_powers = nanstd(rand_theta_powers, 0, 2) ./ sqrt(sum(~isnan(rand_theta_powers), 2));

% Plotting only if there are valid data points
if any(~isnan(mean_powers))
    window_starts = on_ts(1) + (0:windowinc:(max_windows-1)*windowinc);
    valid_indices = ~isnan(rand_mean_powers);
    window_starts = window_starts(valid_indices);
    rand_mean_powers = rand_mean_powers(valid_indices);
    rand_sem_powers = rand_sem_powers(valid_indices);

    figure;
    errorbar(window_starts, mean_powers, rand_sem_powers, 'b-o', 'DisplayName', 'Average Theta Power');
    xlabel('Time (s) relative to Onset');
    ylabel('Average Theta Power (dB) with SEM');
    title('Averaged Theta Power over Time with SEM');
    legend show;
else
    disp('No valid power data was calculated.');
end

% Normalize mean_powers to the highest average theta power value
max_rand_mean_power = max(rand_mean_powers);  % Find the maximum value in mean_powers
normalized_rand_mean_powers = rand_mean_powers / max_rand_mean_power;  % Normalize all values relative to the maximum

