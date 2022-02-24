function tm(parameters)

file = load(parameters.matFileWithSamples);
signal = file.samples;
signal = permute(signal, [2 1 3]);

signal = signal(1:parameters.numberOfAntennas,parameters.tmSamples(1):parameters.tmSamples(2),:);

if parameters.tmBackwardSmoothingUsed
    signal = [signal conj(signal(size(signal,1):-1:1,:,:))];
end

antennasPositions = zeros(parameters.numberOfAntennas - parameters.tmSubarrayNum + 1,3);
for i=1:parameters.numberOfAntennas - parameters.tmSubarrayNum + 1
    antennasPositions(i,1)=parameters.separationBetweenAntennas*i;
end

azimuthAngle = -90:0.1:90; % vector of angle of arrivals with 0.5 deg resolution
azimuthAngle = azimuthAngle.'; % transpose

musicSpectrum = zeros(length(azimuthAngle), parameters.tmNumberOfIterations);
strongestPeaks = zeros(length(azimuthAngle), parameters.tmNumberOfIterations);

for it = 1:parameters.tmNumberOfIterations
    musicSpectrum(:,it) = reallocation(antennasPositions, azimuthAngle, signal(:,:,it), parameters);
    [angles, strongestPeaks(:,it)] = findPeaksOnSpectrum(musicSpectrum(:,it), azimuthAngle, parameters);
    tmEstimatedAngles{it} = angles(angles < 89.5 & angles > -89.5);
end

tmStableStdMean = findstable(tmEstimatedAngles);

figure; plot(azimuthAngle,musicSpectrum), xlabel('Angle, degrees'), ylabel('Amplitude');
title("Time domain MUSIC. Estimated angle: " + string(tmStableStdMean(1,2)));
xlim([-100 100]), ylim([min(musicSpectrum, [], 'all') max(musicSpectrum, [], 'all')]), xticks(-180:45:180);

end % tm

function musicSpectrum = reallocation(antennasPositions, azimuthAngle, signal, parameters)

L = size(antennasPositions,1);

Rxx = zeros(L,L, parameters.tmSubarrayNum);
for i = 1:parameters.tmSubarrayNum
    s = signal(i:i+L-1,:);
    Rxx(:,:,i)=s*s'/(parameters.tmSamples(2)-parameters.tmSamples(1)); % correlation matrix
end
Rxx = sum(Rxx,3) / parameters.tmSubarrayNum;

if parameters.tmMode == 0 % old way, comparison with threshold
    [Utmp, eigenValue]=eig(Rxx);
    eigen=diag(abs(eigenValue));
    [~,I] = sort(eigen, 'ascend');
    U = Utmp(:,I);
    eigen = eigen(I);
    
    N = 0;
    for i=1:size(antennasPositions,1)
        if eigen(i) < parameters.tmThreshold
            N = i;
        end
    end
    
    if N == 0
        N = 1;
    end
    
    En = U(:,1:N);
    
    if N == 0
        error('tm: sources not found')
    end
elseif parameters.tmMode == 1
    [Utmp,eigenValue] = eig(Rxx);
    eigen=diag(abs(eigenValue));
    [~,I] = sort(eigen, 'ascend');
    U = Utmp(:,I);
    eigen = eigen(I);
    
    N = parameters.tmNoiseEigenNumber;
    if N < 0
        N = length(eigen)+N;
    end
    En = U(:,1:N);
end


elevationAngle = zeros(size(azimuthAngle)); % vector of zeros with the same size as az
steeringVectors = steeringVectorsCreate(antennasPositions, [azimuthAngle + 90,elevationAngle], parameters.centralFrequency);
musicSpectrum=zeros(1,length(azimuthAngle));
for i=1:length(azimuthAngle)
    musicSpectrum(i)=-10*log10(abs(steeringVectors(:,i)'*(En*En')*steeringVectors(:,i)));
end
musicSpectrum = musicSpectrum.';

musicSpectrum = musicSpectrum - min(musicSpectrum);

end % reallocation

function steeringVectors = steeringVectorsCreate(antennasPositions, directionOfArrival, centralFrequency)

directionOfArrival=directionOfArrival*pi()/180.0; % Degree to radian conversion
azimuthAngle = directionOfArrival(:,1);
elevationAngle = directionOfArrival(:,2);
projection = [cos(azimuthAngle).*cos(elevationAngle)   sin(azimuthAngle).*cos(elevationAngle)   sin(elevationAngle)]';

waveNumberProjection = 2 * pi() * centralFrequency / physconst('LightSpeed') * projection;
steeringVectors = exp(-1j*(antennasPositions*waveNumberProjection));

end % steeringVectorsCreate

function [tmEstimatedAngles, strongestPeaks] = findPeaksOnSpectrum(musicSpectrum, azimuthAngle, parameters)

switch parameters.tmPeaksMode
    case 0
        [globalMaximumValue, indexOfGlobalMaximum] = max(musicSpectrum);
        tmEstimatedAngles = azimuthAngle(1) + (indexOfGlobalMaximum-1)*(azimuthAngle(2)-azimuthAngle(1));
        strongestPeaks = musicSpectrum .* (musicSpectrum == globalMaximumValue);
    case 1
        isPeak = imregionalmax(musicSpectrum);
        peakIndex = ind2sub(size(isPeak), find(isPeak));
        isStrongEnough = musicSpectrum(isPeak) > (parameters.tmPeakThresholdRate * max(musicSpectrum(isPeak)));
        peakIndex = peakIndex(isStrongEnough);
        tmEstimatedAngles = azimuthAngle(1) + (peakIndex-1)*(azimuthAngle(2)-azimuthAngle(1));
        
        globalMaximumValue = max(musicSpectrum(isPeak));
        pointIsPeak = musicSpectrum & isPeak;
        eachPeakValue = pointIsPeak .*  musicSpectrum;
        isPeakMoreThenThreshold = eachPeakValue > (parameters.tmPeakThresholdRate * globalMaximumValue);
        strongestPeaks = isPeakMoreThenThreshold .* musicSpectrum;
    case 2
        strongestPeaks = zeros(size(musicSpectrum));
        isPeak = imregionalmax(musicSpectrum);
        peakIndex = ind2sub(size(isPeak), find(isPeak));
        [~,topPeakIndices] = maxk(musicSpectrum(isPeak), min(parameters.tmNumberOfPeaksToDetect, length(find(isPeak(:)))));
        peakIndex = peakIndex(topPeakIndices);
        tmEstimatedAngles = azimuthAngle(1) + (peakIndex-1)*(azimuthAngle(2)-azimuthAngle(1));
        
        for i = 1:length(peakIndex)
            strongestPeaks(peakIndex(i)) = musicSpectrum(peakIndex(i));
        end
end
end % findPeaksOnSpectrum