function fm(parameters)

% read CSI from File
subCarrInd = parameters.fmSubCarrIndStart:parameters.fmSubCarrIndStep:parameters.fmSubCarrIndEnd;

file = load(parameters.matFileWithCSI);
matrixCSI = file.csi;
matrixCSI = matrixCSI(subCarrInd, 1:parameters.numberOfAntennas,1:parameters.fmNumberOfIterations * parameters.fmNumberOfPacketsPerIteration);

for mainPartOfComputingSpotfi = 1
    sanitizedVectorCSI = sanitizeCSI(matrixCSI, subCarrInd, parameters);

    sanitizedVectorCSI = permute(sanitizedVectorCSI, [2 1 3]);
    
    if parameters.fmBackwardSmoothingUsed
        sanitizedVectorCSI = [sanitizedVectorCSI conj(sanitizedVectorCSI(size(sanitizedVectorCSI,1):-1:1,:,:))];
    end
    
    L = size(sanitizedVectorCSI, 1) - parameters.fmSubarrayNum + 1;
    Rxx = zeros(L, L, parameters.fmNumberOfIterations);
    for j = 1:parameters.fmNumberOfIterations
        for i = 1:parameters.fmNumberOfPacketsPerIteration
            csi = sanitizedVectorCSI(:,:,i+(j-1)*parameters.fmNumberOfPacketsPerIteration);
            
            for f = 1:parameters.fmSubarrayNum
                s = csi(f:f+L-1,:);
                rxx(:,:,f)=s*s';
            end
            rxx = sum(rxx,3) / parameters.fmSubarrayNum;
    
            Rxx(:,:,j) = Rxx(:,:,j) + rxx; % correlation matrix
        end
    end
    
    Rxx = Rxx / parameters.fmNumberOfPacketsPerIteration;
    
    for j = 1:parameters.fmNumberOfIterations
        if parameters.fmEigenMode == 0 % old way, comparison with threshold
            [Ufmp, eigenValue]=eig(Rxx(:,:,j));
            eigenValue=diag(eigenValue);
            ne=[]; % for noise eigenvectors
            for i=1:length(eigenValue)
                if real(eigenValue(i)) < parameters.fmThreshold
                    ne=[ne i]; % choose vector corresponding to the weakest eigenvalues
                end
            end
            
            if ~isempty(ne)
                En{j}=Ufmp(:,ne);
            else
                numberOfSpatialSources = 1;
                [~,I] = sort(eigenValue, 'ascend');
                U = Ufmp(:,I);
                En{j} = U(:,1:numberOfSpatialSources);
            end
            
        elseif parameters.fmEigenMode == 1
            [Ufmp,eigenValue] = eig(Rxx(:,:,j));
            eigenValue = diag(eigenValue);
            [~,I] = sort(eigenValue, 'ascend');
            U = Ufmp(:,I);
            
            N = parameters.fmNoiseEigenNumber;
            if N < 0
                N = length(eigenValue)+N;
            end
            En{j} = U(:,1:N);
        end
    end
    
    antennasPositions = zeros(parameters.numberOfAntennas - parameters.fmSubarrayNum + 1,3);
    for i=1:size(antennasPositions,1)
        antennasPositions(i,1)=parameters.separationBetweenAntennas*i;
    end
    
    azimuthAngle = [-90:0.1:90] + 90; % vector of angle of arrivals with 0.5 deg resolution
    azimuthAngle = azimuthAngle'; % transpose
    elevationAngle = zeros(size(azimuthAngle)); % vector of zeros with the same size as az
    
    steeringVectors = steeringVectorsCreate(antennasPositions, [azimuthAngle,elevationAngle], parameters.centralFrequency);
    musicSpectrum = zeros(length(azimuthAngle), parameters.fmNumberOfIterations);
    strongestPeaks = zeros(length(azimuthAngle), parameters.fmNumberOfIterations);

    for j = 1:parameters.fmNumberOfIterations
        for i=1:length(azimuthAngle)
            musicSpectrum(i,j)=-10*log10(abs(steeringVectors(:,i)'*(En{j}*En{j}')*steeringVectors(:,i)));
        end
        musicSpectrum(:,j) = musicSpectrum(:,j) - min(musicSpectrum(:,j));
        [angles, strongestPeaks(:,j)] = findPeaksOnSpectrum(musicSpectrum(:,j), azimuthAngle, parameters);
        angles = angles - 90;
        fmEstimatedAngles{j} = angles(angles < 89.5 & angles > -89.5);
    end
    azimuthAngle = azimuthAngle-90;
    
    fmStableStdMean = findstable(fmEstimatedAngles);

    figure; plot(azimuthAngle,musicSpectrum), xlabel('Angle, degrees'), ylabel('Amplitude');
    xlim([-100 100]), ylim([min(musicSpectrum, [], 'all')-2 max(musicSpectrum, [], 'all')+3]), xticks(-180:45:180);
    title("Frequency domain MUSIC. Estimated angle: " + string(fmStableStdMean(1,2)));
end

end

function sanitizedVectorCSI = sanitizeCSI(vectorCSI, subCarrInd, parameters)

sanitizedVectorCSI = zeros(size(vectorCSI));
numberOfAntennas = parameters.numberOfAntennas;
for packetIteration=1:parameters.fmNumberOfIterations*parameters.fmNumberOfPacketsPerIteration
    reshapedCSI = reshape(vectorCSI(:,:,packetIteration), length(subCarrInd), numberOfAntennas);
    
    [PhaseSlope, PhaseConstant] = removePhsSlope(reshapedCSI, subCarrInd, parameters);
    sanitizedVectorCSI(:,:,packetIteration) = reshapedCSI;
end

end

function [PhaseSlope, PhaseConstant] = removePhsSlope(reshapedCSI, subCarrInd, parameters)

PhaseRelatedToFirstPacket = unwrap(angle(reshapedCSI));

for antIdForPhs = 1:parameters.numberOfAntennas
    if  PhaseRelatedToFirstPacket(1,antIdForPhs) - PhaseRelatedToFirstPacket(1,1) > pi
        PhaseRelatedToFirstPacket(:,antIdForPhs) = PhaseRelatedToFirstPacket(:,antIdForPhs) - 2*pi;
    elseif PhaseRelatedToFirstPacket(1,antIdForPhs) - PhaseRelatedToFirstPacket(1,1) < -pi
        PhaseRelatedToFirstPacket(:,antIdForPhs) = PhaseRelatedToFirstPacket(:,antIdForPhs) + 2*pi;
    end
end
A = [repmat(subCarrInd(:), parameters.numberOfAntennas, 1) ones(length(subCarrInd)*parameters.numberOfAntennas, 1)];
x = A\PhaseRelatedToFirstPacket(:);
PhaseSlope = x(1);
PhaseConstant = x(2);

end % removePhsSlope

function steeringVectors = steeringVectorsCreate(antennasPositions, directionOfArrival, centralFrequency)

directionOfArrival=directionOfArrival*pi()/180.0; % Degree to radian conversion
azimuthAngle = directionOfArrival(:,1);
elevationAngle = directionOfArrival(:,2);
projection = [cos(azimuthAngle).*cos(elevationAngle)   sin(azimuthAngle).*cos(elevationAngle)   sin(elevationAngle)]';

waveNumberProjection = 2 * pi() * centralFrequency / physconst('LightSpeed') * projection;
steeringVectors = exp(-1j*(antennasPositions*waveNumberProjection));

end % steeringVectorsCreate

function [EstimatedAngles, strongestPeaks] = findPeaksOnSpectrum(musicSpectrum, azimuthAngle, parameters)

switch parameters.fmPeaksMode
    case 0
        [globalMaximumValue, indexOfGlobalMaximum] = max(musicSpectrum);
        EstimatedAngles = azimuthAngle(1) + (indexOfGlobalMaximum-1)*(azimuthAngle(2)-azimuthAngle(1));
        strongestPeaks = musicSpectrum .* (musicSpectrum == globalMaximumValue);
    case 1
        isPeak = imregionalmax(musicSpectrum);
        peakIndex = ind2sub(size(isPeak), find(isPeak));
        isStrongEnough = musicSpectrum(isPeak) > (parameters.fmPeakThresholdRate * max(musicSpectrum(isPeak)));
        peakIndex = peakIndex(isStrongEnough);
        EstimatedAngles = azimuthAngle(1) + (peakIndex-1)*(azimuthAngle(2)-azimuthAngle(1));
        
        globalMaximumValue = max(musicSpectrum(isPeak));
        pointIsPeak = musicSpectrum & isPeak;
        eachPeakValue = pointIsPeak .*  musicSpectrum;
        isPeakMoreThenThreshold = eachPeakValue > (parameters.fmPeakThresholdRate * globalMaximumValue);
        strongestPeaks = isPeakMoreThenThreshold .* musicSpectrum;
    case 2
        strongestPeaks = zeros(size(musicSpectrum));
        isPeak = imregionalmax(musicSpectrum);
        peakIndex = ind2sub(size(isPeak), find(isPeak));
        [~,topPeakIndices] = maxk(musicSpectrum(isPeak), min(parameters.fmNumberOfPeaksToDetect, length(find(isPeak(:)))));
        peakIndex = peakIndex(topPeakIndices);
        EstimatedAngles = azimuthAngle(1) + (peakIndex-1)*(azimuthAngle(2)-azimuthAngle(1));
        
        for i = 1:length(peakIndex)
            strongestPeaks(peakIndex(i)) = musicSpectrum(peakIndex(i));
        end
end
end % findPeaksOnSpectrum