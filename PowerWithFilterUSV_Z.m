%% Top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior & plot power is USV band
% NOTE: requires Statistics toolbox

%% Initialization

tic                          % track time to run script
clear all; close all; clc;   % reset everything

% set the file path
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\BNST ChR2 USV Master dataset\Male T1 T2 n=19\50Hz\try\';
wavPath = '/Users/zacharygutierrez/Dropbox (Scripps Research)/BNST ChR2 USV Master dataset/Male T1 T2 n=14/try/';

% create a structure of .wav files
files = dir(fullfile(wavPath, '*.wav'));  % fullfile connects the path with .wav   
totalMice = size(files, 1);               % totalMice equals the number of files

AUCpowerARRAY = [];     
nameList = [];
%fig = figure;       % create figure and hold before the loop
%hold on;

usvPowerPerSampleSmoothAll = zeros(1,22463);  % hardcode with amount of validSamples

%% Begin loop for every animal

%totalMice = 5;         % note: temporary toggle for totalMice
for j = 1:totalMice                     
   
waveFile = [wavPath, files(j).name];
disp(['Reading WAV file: ', files(j).name])   % note: can use waveFile instead of files(j).name
nameList = [nameList; {files(j).name(1:end-4)}]; % create a list of animal names without .wav extension
[y, Fs] = audioread(waveFile);

% note: error checking for different sample rates
disp(['Checking the sampling rate: ', num2str(Fs)])
if Fs == 250000
    y = resample(y, 192000, 250000);
end

% [y, Fs]= audioread(files(j).name);
i = audioinfo(waveFile);
totalSec = i.Duration;
startInd = 1;
stopInd = length(y);
yTrim = y(startInd:stopInd);

disp([strcat('Computing FFT #', int2str(j), '...')])
nfft = 512;
window = 512;
noverlap = window*0.5;
% thresh = -90; % threshold in decibels;
[~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh);  %do not use threshold if zscoring below
% note that P is the power spectral density in W/Hz
Tperiod = i.Duration/length(T);

refPower = 10^-12; %reference power is 10?12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 

disp(['Filtering noise...'])
% idea here take z-score & compare to other freqs to remove broadband noise: 
zsignal = zscore(signal);
lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
highFreq = find(F>80000,1,'first'); % index for high cutoff used below
zsignal(1:lowFreq,:) = 0;           % lowpass - set everything below cutoff to zero
zsignal(highFreq:end,:) = 0;        % highpass - add by Jingyi on 11/29/2018
zthresh = 1.5;                      % 1.3 for pup USV
zsignal(zsignal<zthresh) = 0;       % threshold zscore
signalCleaned = signal;             % create a copy to clean below
signalCleaned(zsignal==0) = 0;      % JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!

% calculate power in the whistle snippets over time
disp(['Calculating acoustic power...'])
usvPowerPerSample = mean(signalCleaned); % average power across frequencies (so mean dB from ~40-80kHz, over temporal smoothing filter); do NOT normalize for now

% for filtering/smooothing:
wndsz = round(0.05/Tperiod); % convert seconds to samples
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); % smoothing filter 
AUCpower = trapz(usvPowerPerSample)/totalSec; 
AUCpowerARRAY = [AUCpowerARRAY AUCpower]; 

disp(['Plotting...'])

fig=figure;             %note: had to comment out to make my it work

ax1 = subplot(3,1,1);
%hold(ax1,'on');              % note: add the holds for each axis...
title('original');
imagesc(T, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
set(gca,'YDir','normal')     % flip to make small freq on bottom
colormap(1-bone);
ylabel('frequency (kHz)');
xlabel('time (seconds)')
%hold(ax1,'off');             % note: ...and turn the holds off, repeat

ax2 = subplot(3,1,2);
%hold(ax2,'on');
title('filtered, z-scored, and thresholded');
imagesc(T, F(lowFreq:highFreq)./1000, signalCleaned(lowFreq:highFreq,:)); 
set(gca,'YDir','normal')     % flip to make small freq on bottom
colormap(1-bone);
ylabel('frequency (kHz)');
xlabel('time (seconds)')
%hold(ax2,'off');

ax3 = subplot(3,1,3);
%hold(ax3,'on');
title('power per sample');
plot(ax3, T(1:validSamples), usvPowerPerSampleSmooth)
ylabel({'USV power';'(total dB in 40-80kHz band)'});
xlabel('time (seconds)')
%hold(ax3,'off');

usvPowerPerSampleSmoothAll = cat(1, usvPowerPerSampleSmoothAll,...
     usvPowerPerSampleSmooth);

fprintf(1, '\n');               % put a line break between each sample
% saveas (fig,'tif');
end
%hold off;

AUCAvg = mean(AUCpowerARRAY);

%% After looping and creating averages, plot them

disp(['Loop finished. Now plotting average USV power...']);

% delete first row of zeros or the average suffers
usvPowerPerSampleSmoothAll = usvPowerPerSampleSmoothAll(2:end,:); 

% take the average along the columns
usvPowerPerSampleSmoothAvg = mean(usvPowerPerSampleSmoothAll);

% apply another round of filtering/smoothing
wndsz = round(0.05/Tperiod); % convert seconds to samples;  note: original was 0.05
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
disp((validSamples));
usvPowerPerSampleSmoothAvg2 = conv(usvPowerPerSampleSmoothAvg, gaussFilter, 'same');

figure;
hold on;

% plot the average of the smoothed out signals
plot(T(1:validSamples), usvPowerPerSampleSmoothAvg, 'r:');
title(''); % note: next one overlaps it anyway
ylabel({'USV power';'(total dB in 40-80kHz band)'});
xlabel('time (seconds)')

% plot the smoothed out average (even smoother)
plot(T(1:validSamples), usvPowerPerSampleSmoothAvg2, 'k');
title('Working title...');
ylabel({'USV power';'(total dB in 40-80kHz band)'});
xlabel('time (seconds)')

%plot the stimulus box
xLims = get(gca, 'XLim');
axis([xLims(1) xLims(2) -0.05 2]);    % -0.05 to go slightly below lowest call, 2 should be above highest call
x = [5 10 10 5];
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(1) yLims(2) yLims(2)];
patch(x, y, [0 0.3 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % blue box
set(gca, 'FontSize', 20);     % font size of 20

hold off;

%% Plot the heatmap for Bilateral/L/R/Unilateral
%{
figure;
xLims = get(gca, 'XLim');
colormap('default')
imagesc(usvPowerPerSampleSmoothAll)
title('Sample: EC19');
ylabel({'Trials per stimulus side';'(Bilateral, left, & right; 25/50Hz)'});
xlabel('time (seconds)');
c = colorbar;
c.Label.String = {'USV power';'(total dB in 40-80kHz band)'};
yticklabels({'B25','B25','L25','L25','L25','R25','R25','R25',...
    'B50','B50','L50','L50','L50','R50','R50','R50'});
set(gca, 'YTick', [1:1:16], 'YTickLabel', yticklabels)
xticklabels({'0','5','10','15','20'});
set(gca, 'XTick', [1:4455:22275], 'XTickLabel', xticklabels)
%}

%% Plot the heatmap for Sex/Stimulation frequency

figure;
colormap(flipud(gray))               % flipud swaps the heatmap colors
imagesc(usvPowerPerSampleSmoothAll);
title('Male: 50Hz');
ylabel({'Sample'});
xlabel('time (seconds)');
c = colorbar;
c.Label.String = {'USV power';'(total dB in 40-80kHz band)'};
yticklabels(nameList);
% need as many YTick positions as there are elements in nameList...
set(gca, 'YTick', [1:1:length(nameList)], 'YTickLabel', yticklabels)
set(gca,'TickLabelInterpreter','none');
xticklabels({'0','5','10','15','20','25'});
set(gca, 'XTick', [1:3712.5:22275], 'XTickLabel', xticklabels)

%%
toc % elapsed time to run script              