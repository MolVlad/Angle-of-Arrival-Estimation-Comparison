function parameters = define_parameters

parameters.packetsNumInFile = 30; % for csi.txt and samples.txt
parameters.samplesNumInFile = 1000;
parameters.csiNumInFile = 242;

% files after preprocessing
parameters.matFileWithSamples = "data/samples.mat";
parameters.matFileWithCSI = "data/csi.mat";

parameters.centralFrequency = 2.512e+09; % previously 5.63e9 and 4.997e9
parameters.frequencyGap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
parameters.separationBetweenAntennas = physconst('LightSpeed')/parameters.centralFrequency/2; % lambda/2 - distance between antennas
parameters.numberOfAntennas = 4;   % max 10, after reboot - 14 or 19

parameters.tmNumberOfIterations = 10;
parameters.tmNumberOfPacketsPerIteration = 1;
parameters.tmSamples = [1 1000];

parameters.tmBackwardSmoothingUsed = 0;
parameters.tmSubarrayNum = 1;
parameters.tmMode = 1; % 0 - default with threshold, 1 - fixed amount of noise eig value
parameters.tmPrintEigen = 0;
parameters.tmNoiseEigenNumber = -1; % if < 0 then num of signal
parameters.tmThreshold = 0.1;  % default = 0.001
parameters.tmPeaksMode = 1;  % 0 - absolute max, 1 - threshold from max, 2 - k max peaks
parameters.tmPeakThresholdRate = 0.1; % for tmPeaksMode = 1
parameters.tmNumberOfPeaksToDetect = 8; % for tmPeaksMode = 2
parameters.tmMeasurement = 0;  % to specify is it necessary to write results and dump

parameters.fmNumberOfIterations = 10;
parameters.fmNumberOfPacketsPerIteration = 1;
parameters.fmSubCarrIndStart = 1;
parameters.fmSubCarrIndStep = 1;
parameters.fmSubCarrIndEnd = 242;

parameters.fmBackwardSmoothingUsed = 0;
parameters.fmSubarrayNum = 1;
parameters.fmEigenMode = 1;
parameters.fmThreshold = 0.1;
parameters.fmPeaksMode = 2;  % 0 - absolute max, 1 - threshold from max, 2 - k max peaks
parameters.fmNoiseEigenNumber = -1; % if < 0 then num of signal
parameters.fmPeakThresholdRate = 0.1; % for fmPeaksMode = 1
parameters.fmNumberOfPeaksToDetect = 8; % for fmPeaksMode = 2

end % CreateParams