%% Script to compute local geometric and dose-weighted metrics

clear

%% Loop through all patients
for selectedPatient = 1:250
    clearvars -except selectedPatient metricDosimetric metricGeometric

    %% Load patient data

    % Load patient folders
    patientFolder = ['..\Anonimised data Julian\PAT', num2str(selectedPatient), '\'];
    basisFolder = [patientFolder, 'Basisplan\'];
    fractionFolder = [patientFolder, 'Fr1\'];

    % Load the patient's baseline MR info and data
    [basisMRData, basisMRInfo] = dicomreadVolume([basisFolder, 'MR\']);
    basisMRData = squeeze(basisMRData);

    % Load the patient's baseline delineations
    basisRTstructFilename = dir([basisFolder, 'RTSTRUCT*']).name;
    basisContourData = dicomContours(dicominfo([basisFolder, basisRTstructFilename], UseVRHeuristic = false));

    % Load the patient's baseline dose info and data
    basisRTdoseFilename = dir([basisFolder, 'RTDOSE*']).name;
    basisDoseInfo = dicominfo([basisFolder, basisRTdoseFilename], UseVRHeuristic = false);
    basisDoseDataNormalized = squeeze(dicomread([basisFolder, basisRTdoseFilename], UseVRHeuristic = false));

    basisDoseScalingFactor = basisDoseInfo.DoseGridScaling;
    basisDoseDataAbsolute = basisDoseScalingFactor*double(basisDoseDataNormalized);
    prescribedDose = 36.25;
    basisDoseDataPercentage = basisDoseDataAbsolute/prescribedDose*100;

    % Load the patient's fraction MR info and data
    [fractionMRData, fractionMRInfo] = dicomreadVolume([fractionFolder, 'MR\']);
    fractionMRData = squeeze(fractionMRData);

    % Load the patient's fraction delineations
    fractionRTstructFilename = dir([fractionFolder, 'RTSTRUCT*']).name;
    fractionContourData = dicomContours(dicominfo([fractionFolder, fractionRTstructFilename], UseVRHeuristic = false));

    % Load the patient's auto-fraction-delineations
    autoContourData = importdata([fractionFolder, 'PAT', num2str(selectedPatient), '_predicted_labels.txt']);

    %% Extract baseline image properties

    % Number of rows/columns/slices
    numberOfBasisRows = basisMRInfo.ImageSize(1);
    numberOfBasisColumns = basisMRInfo.ImageSize(2);
    numberOfBasisSlices = basisMRInfo.ImageSize(3);

    % Voxel dimensions
    basisPixelsize = basisMRInfo.PixelSpacings(1);
    basisSlicethickness = basisMRInfo.PatientPositions(2, 3) - basisMRInfo.PatientPositions(1, 3);

    % Coordinates of center of voxel in first row/column/slice
    basisFirstvoxelCenterX = basisMRInfo.PatientPositions(1, 1);
    basisFirstvoxelCenterY = basisMRInfo.PatientPositions(1, 2);
    basisFirstvoxelCenterZ = basisMRInfo.PatientPositions(1, 3);

    % Baseline MR image range
    basisXMin = basisFirstvoxelCenterX - basisPixelsize/2;
    basisXMax = basisXMin + numberOfBasisColumns*basisPixelsize;
    basisYMin = basisFirstvoxelCenterY - basisPixelsize/2;
    basisYMax = basisYMin + numberOfBasisRows*basisPixelsize;
    basisZMin = basisFirstvoxelCenterZ - basisSlicethickness/2;
    basisZMax = basisZMin + numberOfBasisSlices*basisSlicethickness;

    % Isocenter location in baseline image
    basisIsocROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'Isocenter');
    basisIsocCoordinates = cell2mat(basisContourData.ROIs{basisIsocROIRownumber, 3}{1, 1});
    basisIsocVoxel = ...
        [0.5 + (basisIsocCoordinates(2) - basisYMin)/basisPixelsize, ...
        0.5 + (basisIsocCoordinates(1) - basisXMin)/basisPixelsize, ...
        0.5 + (basisIsocCoordinates(3) - basisZMin)/basisSlicethickness];

    %% Extract fraction image properties

    numberOfFractionRows = fractionMRInfo.ImageSize(1);
    numberOfFractionColumns = fractionMRInfo.ImageSize(2);
    numberOfFractionSlices = fractionMRInfo.ImageSize(3);

    fractionPixelsize = fractionMRInfo.PixelSpacings(1);
    fractionSlicethickness = fractionMRInfo.PatientPositions(2, 3) - fractionMRInfo.PatientPositions(1, 3);

    fractionFirstvoxelCenterX = fractionMRInfo.PatientPositions(1, 1);
    fractionFirstvoxelCenterY = fractionMRInfo.PatientPositions(1, 2);
    fractionFirstvoxelCenterZ = fractionMRInfo.PatientPositions(1, 3);

    fractionXMin = fractionFirstvoxelCenterX - fractionPixelsize/2;
    fractionXMax = fractionXMin + numberOfFractionColumns*fractionPixelsize;
    fractionYMin = fractionFirstvoxelCenterY - fractionPixelsize/2;
    fractionYMax = fractionYMin + numberOfFractionRows*fractionPixelsize;
    fractionZMin = fractionFirstvoxelCenterZ - fractionSlicethickness/2;
    fractionZMax = fractionZMin + numberOfFractionSlices*fractionSlicethickness;

    fractionIsocROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'Isocenter');
    fractionIsocCoordinates = cell2mat(fractionContourData.ROIs{fractionIsocROIRownumber, 3}{1, 1});
    fractionIsocVoxel = ...
        [0.5 + (fractionIsocCoordinates(2) - fractionYMin)/fractionPixelsize, ...
        0.5 + (fractionIsocCoordinates(1) - fractionXMin)/fractionPixelsize, ...
        0.5 + (fractionIsocCoordinates(3) - fractionZMin)/fractionSlicethickness];

    %% Match baseline CTV and dose distribution onto fraction image

    % Shift in x,y,z coordinates and in number of rows/columns/slices
    shiftCoordinates = fractionIsocCoordinates - basisIsocCoordinates;
    shiftVoxels = fractionIsocVoxel - basisIsocVoxel;

    % Shifting baseline CTV contour point coordinates
    basisCTVROIRownumber = find(string(basisContourData.ROIs{:, 2}) == 'CTV');
    basisCTVContourpoints = basisContourData.ROIs{basisCTVROIRownumber, 3}{1, 1};
    numberOfBasisCTVContours = length(basisCTVContourpoints);
    for counter_basiscontours = 1:numberOfBasisCTVContours
        basisCTVContourpointsShifted{counter_basiscontours, 1} = ...
            basisCTVContourpoints{counter_basiscontours, 1} + shiftCoordinates;
    end

    % Shifting baseline dose distribution voxel numbers
    basisRows = 1:numberOfBasisRows;
    basisColumns = 1:numberOfBasisColumns;
    basisSlices = 1:numberOfBasisSlices;
    [basisColumnNumbers, basisRowNumbers, basisSliceNumbers] = ...
        meshgrid(basisColumns, basisRows, basisSlices);

    matchingFractionRowNumbers = basisRowNumbers + shiftVoxels(1);
    matchingFractionColumnNumbers = basisColumnNumbers + shiftVoxels(2);
    matchingFractionSliceNumbers = basisSliceNumbers + shiftVoxels(3);

    basisDoseDataAbsoluteShifted = ...
        interp3(matchingFractionColumnNumbers, matchingFractionRowNumbers, matchingFractionSliceNumbers, ...
        basisDoseDataAbsolute, basisColumnNumbers, basisRowNumbers, basisSliceNumbers);

    %% Save matched baseline CTV as ROI

    matchedCTVContourData = deleteContour(basisContourData, 1:basisIsocROIRownumber);
    matchedCTVContourData = addContour(matchedCTVContourData, 1, ...
        'CTV_basis_matched', basisCTVContourpointsShifted, 'Closed_planar');

    %% Pixelize baseline- and auto-fraction CTV delineations on fraction images
    
    fractionMRSpatialInfo = imref3d([numberOfFractionRows, numberOfFractionColumns, numberOfFractionSlices], ...
        [fractionXMin + fractionPixelsize/2, fractionXMax + fractionPixelsize/2], ...
        [fractionYMin + fractionPixelsize/2, fractionYMax + fractionPixelsize/2], ...
        [fractionZMin + fractionSlicethickness/2, fractionZMax + fractionSlicethickness/2]);
    
    % Mask matched baseline CTV delineation
    basisMaskTrans = createMask(matchedCTVContourData, 'CTV_basis_matched', fractionMRSpatialInfo);
    
    % Mask cropped auto-fraction delineations
    autoCroppedRows = 63:238;
    autoCroppedColumns = 32:303;
    autoCroppedSlices = 89:200;
    autoMasksTransCropped = flip(permute(reshape(autoContourData, [272, 176, 112]), [2, 1, 3]), 3);
    autoMasksTransFull = zeros(numberOfFractionRows, numberOfFractionColumns, numberOfFractionSlices);
    autoMasksTransFull(autoCroppedRows, autoCroppedColumns, autoCroppedSlices) = autoMasksTransCropped;
    % Only keep CTV delineation
    autoMaskTrans = autoMasksTransFull;
    autoMaskTrans(autoMaskTrans ~= 1) = 0;

    %% Extract CTV masks in remaining planes for 3D boundary voxel computation

    autoMaskCor = permute(autoMaskTrans, [3, 2, 1]);
    autoMaskSag = permute(autoMaskTrans, [3, 1, 2]);
    basisMaskCor = permute(basisMaskTrans, [3, 2, 1]);
    basisMaskSag = permute(basisMaskTrans, [3, 1, 2]);

    %% Extract boundary voxel numbers from transversal masks

    for counter_slices = 1:numberOfFractionSlices
        % Generate transversal boundary images
        autoTransBoundaryImage(:, :, counter_slices) = bwmorph(autoMaskTrans(:, :, counter_slices), 'remove');
        basisTransBoundaryImage(:, :, counter_slices) = bwmorph(basisMaskTrans(:, :, counter_slices), 'remove');
    
        % Extract transversal boundary pixel numbers (row, column)
        [autoTransBoundaryVoxelRowNumbersBySlice{counter_slices, 1}, ...
            autoTransBoundaryVoxelColumnNumbersBySlice{counter_slices, 1}] = ...
            ind2sub([numberOfFractionRows, numberOfFractionColumns], find(autoTransBoundaryImage(:, :, counter_slices) == 1));
        % Complete with slice numbers
        autoTransBoundaryVoxelSliceNumbersBySlice{counter_slices, 1} = ...
            counter_slices*ones(length(autoTransBoundaryVoxelRowNumbersBySlice{counter_slices, 1}), 1);

        [basisTransBoundaryVoxelRowNumbersBySlice{counter_slices, 1}, ...
            basisTransBoundaryVoxelColumnNumbersBySlice{counter_slices, 1}] = ...
            ind2sub([numberOfFractionRows, numberOfFractionColumns], find(basisTransBoundaryImage(:, :, counter_slices) == 1));
        basisTransBoundaryVoxelSliceNumbersBySlice{counter_slices, 1} = ...
            counter_slices*ones(length(basisTransBoundaryVoxelRowNumbersBySlice{counter_slices, 1}), 1);

        % Convert to voxel numbers
        autoTransBoundaryVoxelNumbersBySlice{counter_slices, 1} = ...
            cat(2, autoTransBoundaryVoxelRowNumbersBySlice{counter_slices, 1}, ...
            autoTransBoundaryVoxelColumnNumbersBySlice{counter_slices, 1}, ...
            autoTransBoundaryVoxelSliceNumbersBySlice{counter_slices, 1}); 

        basisTransBoundaryVoxelNumbersBySlice{counter_slices, 1} = ...
            cat(2, basisTransBoundaryVoxelRowNumbersBySlice{counter_slices, 1}, ...
            basisTransBoundaryVoxelColumnNumbersBySlice{counter_slices, 1}, ...
            basisTransBoundaryVoxelSliceNumbersBySlice{counter_slices, 1}); 
    end

    % Save all boundary voxels numbers in a single cell
    autoTransBoundaryVoxelNumbers = cat(1, autoTransBoundaryVoxelNumbersBySlice{:});
    basisTransBoundaryVoxelNumbers = cat(1, basisTransBoundaryVoxelNumbersBySlice{:});

    %% Same for coronal masks

    for counter_rows = 1:numberOfFractionRows
        autoCorBoundaryImage(:, :, counter_rows) = bwmorph(autoMaskCor(:, :, counter_rows), 'remove');
        basisCorBoundaryImage(:, :, counter_rows) = bwmorph(basisMaskCor(:, :, counter_rows), 'remove');
    
        [autoCorBoundaryVoxelSliceNumbersByRow{counter_rows, 1}, ...
            autoCorBoundaryVoxelColumnNumbersByRow{counter_rows, 1}] = ...
            ind2sub([numberOfFractionSlices, numberOfFractionColumns], find(autoCorBoundaryImage(:, :, counter_rows) == 1));
        autoCorBoundaryVoxelRowNumbersByRow{counter_rows, 1} = ...
            counter_rows*ones(length(autoCorBoundaryVoxelSliceNumbersByRow{counter_rows, 1}), 1);

        [basisCorBoundaryVoxelSliceNumbersByRow{counter_rows, 1}, ...
            basisCorBoundaryVoxelColumnNumbersByRow{counter_rows, 1}] = ...
            ind2sub([numberOfFractionSlices, numberOfFractionColumns], find(basisCorBoundaryImage(:, :, counter_rows) == 1));
        basisCorBoundaryVoxelRowNumbersByRow{counter_rows, 1} = ...
            counter_rows*ones(length(basisCorBoundaryVoxelSliceNumbersByRow{counter_rows, 1}), 1);
    
        autoCorBoundaryVoxelNumbersByRow{counter_rows, 1} = ...
            cat(2, autoCorBoundaryVoxelRowNumbersByRow{counter_rows, 1}, ...
            autoCorBoundaryVoxelColumnNumbersByRow{counter_rows, 1}, ...
            autoCorBoundaryVoxelSliceNumbersByRow{counter_rows, 1});

        basisCorBoundaryVoxelNumbersByRow{counter_rows, 1} = ...
            cat(2, basisCorBoundaryVoxelRowNumbersByRow{counter_rows, 1}, ...
            basisCorBoundaryVoxelColumnNumbersByRow{counter_rows, 1}, ...
            basisCorBoundaryVoxelSliceNumbersByRow{counter_rows, 1}); 
    end

    autoCorBoundaryVoxelNumbers = cat(1, autoCorBoundaryVoxelNumbersByRow{:});
    basisCorBoundaryVoxelNumbers = cat(1, basisCorBoundaryVoxelNumbersByRow{:});

    %% And same for saggital masks

    for counter_columns = 1:numberOfFractionColumns
        autoSagBoundaryImage(:, :, counter_columns) = bwmorph(autoMaskSag(:, :, counter_columns), 'remove');
        basisSagBoundaryImage(:, :, counter_columns) = bwmorph(basisMaskSag(:, :, counter_columns), 'remove');
    
        [autoSagBoundaryVoxelSliceNumbersByColumn{counter_columns, 1}, ...
            autoSagBoundaryVoxelRowNumbersByColumn{counter_columns, 1}] = ...
            ind2sub([numberOfFractionSlices, numberOfFractionColumns], find(autoSagBoundaryImage(:, :, counter_columns) == 1));
        autoSagBoundaryVoxelColumnNumbersByColumn{counter_columns, 1} = ...
            counter_columns*ones(length(autoSagBoundaryVoxelSliceNumbersByColumn{counter_columns, 1}), 1);

        [basisSagBoundaryVoxelSliceNumbersByColumn{counter_columns, 1}, ...
            basisSagBoundaryVoxelRowNumbersByColumn{counter_columns, 1}] = ...
            ind2sub([numberOfFractionSlices, numberOfFractionColumns], find(basisSagBoundaryImage(:, :, counter_columns) == 1));
        basisSagBoundaryVoxelColumnNumbersByColumn{counter_columns, 1} = ...
            counter_columns*ones(length(basisSagBoundaryVoxelSliceNumbersByColumn{counter_columns, 1}), 1);

        autoSagBoundaryVoxelNumbersByColumn{counter_columns, 1} = ...
            cat(2, autoSagBoundaryVoxelRowNumbersByColumn{counter_columns, 1}, ...
            autoSagBoundaryVoxelColumnNumbersByColumn{counter_columns, 1}, ...
            autoSagBoundaryVoxelSliceNumbersByColumn{counter_columns, 1}); 

        basisSagBoundaryVoxelNumbersByColumn{counter_columns, 1} = ...
            cat(2, basisSagBoundaryVoxelRowNumbersByColumn{counter_columns, 1}, ...
            basisSagBoundaryVoxelColumnNumbersByColumn{counter_columns, 1}, ...
            basisSagBoundaryVoxelSliceNumbersByColumn{counter_columns, 1}); 
    end

   autoSagBoundaryVoxelNumbers = cat(1, autoSagBoundaryVoxelNumbersByColumn{:});
    basisSagBoundaryVoxelNumbers = cat(1, basisSagBoundaryVoxelNumbersByColumn{:});

    %% Combine boundary voxels in 2D planes to extract boundary voxels of 3D objects

    autoBoundaryVoxelNumbers = unique(cat(1, ...
    autoCorBoundaryVoxelNumbers, autoSagBoundaryVoxelNumbers, autoTransBoundaryVoxelNumbers), 'rows');
    numberOfAutoBoundaryVoxels = length(autoBoundaryVoxelNumbers);

    basisBoundaryVoxelNumbers = unique(cat(1, ...
        basisCorBoundaryVoxelNumbers, basisSagBoundaryVoxelNumbers, basisTransBoundaryVoxelNumbers), 'rows');
    numberOfBasisBoundaryVoxels = length(basisBoundaryVoxelNumbers);

    %% Convert boundary voxel numbers to Cartesian coordinates

    basisBoundaryVoxelCoordinates = ...
        [(basisBoundaryVoxelNumbers(:,2) - fractionIsocVoxel(2))*fractionPixelsize, ...
        (basisBoundaryVoxelNumbers(:,1) - fractionIsocVoxel(1))*fractionPixelsize, ...
        (basisBoundaryVoxelNumbers(:,3) - fractionIsocVoxel(3))*fractionSlicethickness];

    autoBoundaryVoxelCoordinates = ...
        [(autoBoundaryVoxelNumbers(:,2) - fractionIsocVoxel(2))*fractionPixelsize, ...
        (autoBoundaryVoxelNumbers(:,1) - fractionIsocVoxel(1))*fractionPixelsize, ...
        (autoBoundaryVoxelNumbers(:,3) - fractionIsocVoxel(3))*fractionSlicethickness];

    %% Center coordinates such that origin lies at centroid of baseline CTV delineation

    originCoordinates = [mean(basisBoundaryVoxelCoordinates(:,1)), 
        mean(basisBoundaryVoxelCoordinates(:,2)), mean(basisBoundaryVoxelCoordinates(:,3))];

    basisBoundaryVoxelCoordinatesCentred = ...
        [(basisBoundaryVoxelCoordinates(:,1) - originCoordinates(1)), ...
        (basisBoundaryVoxelCoordinates(:,2) - originCoordinates(2)), ...
        (basisBoundaryVoxelCoordinates(:,3) - originCoordinates(3))];

    autoBoundaryVoxelCoordinatesCentred = ...
        [(autoBoundaryVoxelCoordinates(:,1) - originCoordinates(1)), ...
        (autoBoundaryVoxelCoordinates(:,2) - originCoordinates(2)), ...
        (autoBoundaryVoxelCoordinates(:,3) - originCoordinates(3))];

    %% Convert centred Cartesian coordinates in baseline CTV to spherical coordinates

    % Conversion from Cartesian to Matlab's convention of spherical coordinates
    [azimuthalCoordinate, elevationCoordinate, rCoordinate] = cart2sph(basisBoundaryVoxelCoordinatesCentred(:,1), ...
        basisBoundaryVoxelCoordinatesCentred(:,2), basisBoundaryVoxelCoordinatesCentred(:,3));

    % Conversion from Matlab's convention to standard spherical coordinates
    % Transform phi from (-pi, pi) to (0, 2pi)
    % Transform theta from theta_z=0 = 0 to theta_z=0 = pi/2
    for counter_basisboundaryvoxels = 1:numberOfBasisBoundaryVoxels
        if azimuthalCoordinate(counter_basisboundaryvoxels) < 0
            phiCoordinate(counter_basisboundaryvoxels) = azimuthalCoordinate(counter_basisboundaryvoxels) + 2*pi;
        else phiCoordinate(counter_basisboundaryvoxels) = azimuthalCoordinate(counter_basisboundaryvoxels);
        end
    end

    phiCoordinate = phiCoordinate';
    thetaCoordinate = pi/2 - elevationCoordinate;

    %% Divide basis boundary voxels into segments

    numberOfThetaIntervals = 5;
    numberOfPhiIntervals = 4;
    numberOfSegments = numberOfThetaIntervals*numberOfPhiIntervals;

    % Extract segment-specific indices of basis boundary voxels
    for counter_deltatheta = 1:numberOfThetaIntervals
        for counter_deltaphi = 1:numberOfPhiIntervals
            basisBoundaryVoxelSegmentIndices{counter_deltatheta, counter_deltaphi} = find( ...
                (counter_deltatheta - 1)*pi/numberOfThetaIntervals < thetaCoordinate & ...
                thetaCoordinate < counter_deltatheta*pi/numberOfThetaIntervals & ...
                (counter_deltaphi - 1)*2*pi/numberOfPhiIntervals < phiCoordinate & ...
                phiCoordinate < counter_deltaphi*2*pi/numberOfPhiIntervals);
        end
    end

    % Save segment-specific indices in single column
    basisBoundaryVoxelSegmentIndicesSingleColumn = reshape(basisBoundaryVoxelSegmentIndices', [numberOfSegments, 1]);

    % Extract corresponding basis boundary voxel coordinates 
    for counter_segments = 1:numberOfSegments
        basisBoundaryVoxelCoordinatesSegmented{counter_segments, 1} = ...
            basisBoundaryVoxelCoordinatesCentred(basisBoundaryVoxelSegmentIndicesSingleColumn{counter_segments, 1}, :);
    end

    %% Compute 3D local editing vectors

    % Find Euclidean distance between all baseline-auto CTV boundary voxel pairs
    voxelDistances = pdist2(basisBoundaryVoxelCoordinatesCentred, autoBoundaryVoxelCoordinatesCentred);

    for counter_basisboundaryvoxels2 = 1:numberOfBasisBoundaryVoxels
        % Save startpoint coordinates in cells
        LEVStartpoints{counter_basisboundaryvoxels2, 1} = ...
            basisBoundaryVoxelCoordinatesCentred(counter_basisboundaryvoxels2, :);

        % Find shortest vector(s) from point on baseline boundary to point on auto boundary
        LEVEndpointIndices{counter_basisboundaryvoxels2, 1} = ...
            find(voxelDistances(counter_basisboundaryvoxels2, :) == min(voxelDistances(counter_basisboundaryvoxels2, :)));
    
        % Find number of local LEVs
        numberOfLocalLEVs(counter_basisboundaryvoxels2, 1) = length(LEVEndpointIndices{counter_basisboundaryvoxels2, 1});
    
        % Save endpoint(s) coordinates in cells
        LEVEndpoints{counter_basisboundaryvoxels2, 1} = ...
            autoBoundaryVoxelCoordinatesCentred(LEVEndpointIndices{counter_basisboundaryvoxels2, 1}, :);
    
        % Subtract startpoint from endpoint(s) for vector components
        LEVComponents{counter_basisboundaryvoxels2, 1} = ...
            LEVEndpoints{counter_basisboundaryvoxels2, 1} - LEVStartpoints{counter_basisboundaryvoxels2, 1};
    end

    %% Compute local baseline dose gradients (in Gy/mm), 
    %  take inner product with LEVs, and keep largest local values

    [basisDoseGradientX, basisDoseGradientY, basisDoseGradientZ] = ...
        gradient(basisDoseDataAbsoluteShifted, fractionPixelsize, fractionPixelsize, fractionSlicethickness);

    % Find indices of dose gradient at boundary voxels
    basisBoundaryVoxelIndices = sub2ind([numberOfFractionRows, numberOfFractionColumns, numberOfFractionSlices], ...
        basisBoundaryVoxelNumbers(:, 1), basisBoundaryVoxelNumbers(:, 2), basisBoundaryVoxelNumbers(:, 3));

    % Extract local baseline dose gradients at baseline boundary voxels
    basisDoseGradientXAtBasisBoundary = basisDoseGradientX(basisBoundaryVoxelIndices);
    basisDoseGradientYAtBasisBoundary = basisDoseGradientY(basisBoundaryVoxelIndices);
    basisDoseGradientZAtBasisBoundary = basisDoseGradientZ(basisBoundaryVoxelIndices);
    basisDoseGradientComponentsAtBasisBoundary = [basisDoseGradientXAtBasisBoundary, ...
        basisDoseGradientYAtBasisBoundary, basisDoseGradientZAtBasisBoundary];

    for counter_basisboundaryvoxels3 = 1:numberOfBasisBoundaryVoxels
        % Save dose gradient components in cells
        basisDoseGradientComponentsAtBasisBoundaryCells{counter_basisboundaryvoxels3, 1} = ...
            basisDoseGradientComponentsAtBasisBoundary(counter_basisboundaryvoxels3, :);
    
        % Copy dose gradient components for boundary voxels with multiple LEVs
        basisDoseGradientComponentsAtBasisBoundaryCellsWithDuplicates{counter_basisboundaryvoxels3, 1} = ...
            repmat(basisDoseGradientComponentsAtBasisBoundaryCells{counter_basisboundaryvoxels3, 1}, ...
            numberOfLocalLEVs(counter_basisboundaryvoxels3), 1);

        % Take absolute value of inner product with LEVs
        LEVxDoseGradientAbsolute{counter_basisboundaryvoxels3, 1} = abs(dot(LEVComponents{counter_basisboundaryvoxels3, 1}, ...
            basisDoseGradientComponentsAtBasisBoundaryCellsWithDuplicates{counter_basisboundaryvoxels3, 1}, 2));

        % Find local LEV indices to keep (largest inner product)
        LocalLEVIndicesToKeep{counter_basisboundaryvoxels3, 1} = ...
            find(LEVxDoseGradientAbsolute{counter_basisboundaryvoxels3, 1} == ...
            max(LEVxDoseGradientAbsolute{counter_basisboundaryvoxels3, 1}), 1);

        % Only save selected indices
        LEVComponentsSelected{counter_basisboundaryvoxels3, 1} = ...
            LEVComponents{counter_basisboundaryvoxels3, 1}(LocalLEVIndicesToKeep{counter_basisboundaryvoxels3, 1}, :);

        LEVxDoseGradientAbsoluteSelected{counter_basisboundaryvoxels3, 1} = ...
            LEVxDoseGradientAbsolute{counter_basisboundaryvoxels3, 1}(LocalLEVIndicesToKeep{counter_basisboundaryvoxels3, 1}, :);
    end

    % Save updated variables in single cells
    LEVComponentsSelectedSingleCell = cat(1, LEVComponentsSelected{:});
    LEVxDoseGradientAbsoluteSelectedSingleCell = cat(1, LEVxDoseGradientAbsoluteSelected{:});

    %% Compute segment-specific metrics

    % Input for purely geometric metric
    LEVMagnitude = sqrt(sum(LEVComponentsSelectedSingleCell.^2, 2));

    for counter_segments2 = 1:numberOfSegments
        % Divide inputs over segments
        LEVxDoseGradientAbsoluteSegmented{counter_segments2, 1} = ...
            LEVxDoseGradientAbsoluteSelectedSingleCell(basisBoundaryVoxelSegmentIndicesSingleColumn{counter_segments2, 1});

        LEVMagnitudeSegmented{counter_segments2, 1} = ...
            LEVMagnitude(basisBoundaryVoxelSegmentIndicesSingleColumn{counter_segments2, 1});
    
        % Take mean over all points within segment
        metricDosimetric(selectedPatient, counter_segments2) =  mean(LEVxDoseGradientAbsoluteSegmented{counter_segments2, 1});
        metricGeometric(selectedPatient, counter_segments2) = mean(LEVMagnitudeSegmented{counter_segments2, 1});
    end
end   