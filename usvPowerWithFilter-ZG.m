%% Top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior & plot power is USV band
% NOTE: Requires Statistics Toolbox
% LAST UPDATED: March 12th, 2019 @ 2300

%% Initialization

tic % track time to run script
clear all; close all; clc; % reset everything

% Set the file path for Windows
%wavPath = 'E:\Zach\Dropbox (Scripps Research)\BNST ChR2 USV Master dataset\Male T1 T2 n=19\10Hz\';

% Set the file path for Mac
wavPath = '/Users/zacharygutierrez/Dropbox (Scripps Research)/BNST ChR2 USV Master dataset/Male T1 T2 n=19/1Hz/';

% create a structure of .wav files
files = dir(fullfile(wavPath, '*.wav'));  % fullfile connects the path with .wav   
totalMice = size(files, 1);               % totalMice equals the number of files

AUCpowerARRAY = [];     
nameList = [];
usvPowerPerSampleSmoothMat = zeros(1,22463);  % hardcode with amount of validSamples

%% Begin loop for every animal

%totalMice = 3; % quick toggle for debug

for j = 1:totalMice                     
   
waveFile = [wavPath, files(j).name];
disp(['Reading WAV file: ', files(j).name])   % note: can use waveFile instead of files(j).name
nameList = [nameList; {files(j).name(1:end-4)}]; % create a list of animal names without .wav extension
[y, Fs] = audioread(waveFile);

% error checking for different sample rates -ZG
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
%{
disp(['Plotting...'])

fig=figure;             % note: had to comment out to make mine work

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
%}
usvPowerPerSampleSmoothMat = cat(1, usvPowerPerSampleSmoothMat,...
     usvPowerPerSampleSmooth);

fprintf(1, '\n');               % put a line break between each sample
% saveas (fig,'tif');
end

AUCAvg = mean(AUCpowerARRAY);

%% Plot the average signal -ZG

disp(['Loop finished. Now plotting average USV power...']);

% delete first row of zeros or the average suffers
usvPowerPerSampleSmoothMat = usvPowerPerSampleSmoothMat(2:end,:); 

% take the average along the columns
usvPowerPerSampleSmoothAvg = mean(usvPowerPerSampleSmoothMat);

% apply another round of filtering/smoothing
wndsz = round(0.05/Tperiod);   % convert seconds to samples, original was 0.05
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
usvPowerPerSampleSmoothAvg2 = conv(usvPowerPerSampleSmoothAvg, gaussFilter, 'same');

figure;
hold on;

% plot the average of the smoothed out signals
plot(T(1:validSamples), usvPowerPerSampleSmoothAvg, 'r:');
title(''); % note: next one overlaps it anyway
ylabel({'USV power';'(total dB in 40-80kHz band)'});
xlabel('time (seconds)')

% plot the smoothed out average (even smoother)
% plot(T(1:validSamples), usvPowerPerSampleSmoothAvg2, 'k');
% title('Working title...');
% ylabel({'USV power';'(total dB in 40-80kHz band)'});
% xlabel('time (seconds)')

%plot the stimulus box
xLims = get(gca, 'XLim');
axis([xLims(1) xLims(2) -0.05 2]);    % -0.05 to go slightly below lowest call, 2 should be above highest call
x = [5 10 10 5];
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(1) yLims(2) yLims(2)];
patch(x, y, [0 0.3 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % blue box
set(gca, 'FontSize', 20);     % font size of 20

hold off;

%% Plot the heatmap for Bilateral/L/R/Unilateral -ZG
%{
figure;
xLims = get(gca, 'XLim');
colormap('default')
imagesc(usvPowerPerSampleSmoothMat)
title('Sample: EC19');
ylabel({'Trials per stimulus side';'(Bilateral, left, & right; 25/50Hz)'});
xlabel('time (seconds)');
c = colorbar;
c.Label.String = {'USV power';'(total dB in 40-80kHz band)'};
yticklabels({'B25','B25','L25','L25','L25','R25','R25','R25',...
    'B50','B50','L50','L50','L50','R50','R50','R50'});
set(gca, 'YTick', [1:1:16], 'YTickLabel', yticklabels)
xticklabels({'0','5','10','15','20','25'});
set(gca, 'XTick', [1:3712.5:22275], 'XTickLabel', xticklabels)
%}

%% Sort the matrix rows in descending order of power -ZG

Arowmax = max(usvPowerPerSampleSmoothMat, [], 2);
[~,idx] = sort(Arowmax, 'descend');
SortedMatrix = usvPowerPerSampleSmoothMat(idx,:);
SortedNameList = nameList(idx); % update the index as well
usvPowerPerSampleSmoothMat = SortedMatrix;

%% Plot the heatmap for Sex/Stimulation frequency -ZG

figure;
colormap(flipud(gray))               % flipud swaps the heatmap colors
imagesc(usvPowerPerSampleSmoothMat); % throw a matrix into imagesc, boom heatmap

x=[3712.5 3712.5];   % plot the stimulus start at 5 seconds (3712.5)           
y=[0 totalMice+1];   % should run from bottom (0) to top (totalMice)
line(x,y, 'Color', [0 0.3 1])
x=[7425 7425];       % plot the stimulus end at 10 seconds (7425)           
y=[0 totalMice+1];   % should run from bottom (0) to top (totalMice)
line(x,y, 'Color', [0 0.3 1])

title('Male: 1Hz');
ylabel({'Sample'});
xlabel('time (seconds)');
c = colorbar;
c.Label.String = {'USV power';'(total dB in 40-80kHz band)'};
yticklabels(SortedNameList); 
% need as many YTick positions as there are elements in the name list
set(gca, 'YTick', [1:1:length(SortedNameList)], 'YTickLabel', yticklabels)
set(gca,'TickLabelInterpreter','none');
xticklabels({'0','5','10','15','20','25'});
set(gca, 'XTick', [1:3712.5:22275], 'XTickLabel', xticklabels)

%% Elapsed time to run script -ZG
toc            