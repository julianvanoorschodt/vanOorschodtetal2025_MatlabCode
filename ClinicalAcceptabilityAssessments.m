%% Script to generate GUI that radiation oncologist is provided with for clinical acceptability assessments

clear

% Set default interpreters to Latex (still need to set text interpreters in uilabels to Latex manually)
set(groot,'defaultLegendInterpreter', 'latex');
set(groot,'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex')

% Create empty uifigure on CZE desktop fullscreen size
CZEFullwidth = 1921;
CZEFullheight = 1010;
f1 = uifigure('Position', [0, 42, CZEFullwidth, CZEFullheight]);

% Take a short break to ensure fullscreen figure has popped up
pause(2)

%% Select patient and load patient data

% Ask user to select a patient
inputTitle = 'Patient selection';
inputPrompt = {'Select a patient: (number between 1 and 250)'};
inputValue = inputdlg(inputPrompt, inputTitle, [1 45]);

selectedPatient = str2double(inputValue{1});

% Define valid user input
validPatientValues = 1:250;

% Show error until user input is valid
while ismember(selectedPatient, validPatientValues) == 0
    uialert(f1, 'Please select a valid patient number between 1 and 250!', ...
        'Patient selection error', 'CloseFcn', 'uiresume(f1)');
    uiwait(f1)

    inputValue = inputdlg(inputPrompt, inputTitle, [1 45]);
    selectedPatient = str2double(inputValue{1});
end

% Start loading patient data once input is valid
loadingBar = waitbar(0, 'Searching for patient folder...');
pause(1)

% Load the folders that contain the patient's data
patientFolder = ['..\Anonimised data Julian\PAT', num2str(selectedPatient), '\'];
basisFolder = [patientFolder, 'Basisplan\'];
fractionFolder = [patientFolder, 'Fr1\'];

% Load the patient's baseline MR info and data
waitbar(1/16, loadingBar, 'Loading baseline data...')
[basisMRData, basisMRInfo] = dicomreadVolume([basisFolder, 'MR\']);
basisMRData = squeeze(basisMRData);

% Load the patient's baseline delineations
waitbar(2/16, loadingBar)
basisRTstructFilename = dir([basisFolder, 'RTSTRUCT*']).name;
basisContourData = dicomContours(dicominfo([basisFolder, basisRTstructFilename], UseVRHeuristic = false));

% Load the patient's baseline dose info and data
waitbar(3/16, loadingBar)
basisRTdoseFilename = dir([basisFolder, 'RTDOSE*']).name;
basisDoseInfo = dicominfo([basisFolder, basisRTdoseFilename], UseVRHeuristic = false);
basisDoseDataNormalized = squeeze(dicomread([basisFolder, basisRTdoseFilename], UseVRHeuristic = false));

waitbar(4/16, loadingBar)
basisDoseScalingFactor = basisDoseInfo.DoseGridScaling;
basisDoseDataAbsolute = basisDoseScalingFactor*double(basisDoseDataNormalized);
prescribedDose = 36.25;
basisDoseDataPercentage = basisDoseDataAbsolute/prescribedDose*100;

% Load the patient's fraction MR info and data
waitbar(5/16, loadingBar, 'Loading fraction data...');
[fractionMRData, fractionMRInfo] = dicomreadVolume([fractionFolder, 'MR\']);
fractionMRData = squeeze(fractionMRData);

% Load the patient's fraction delineations
waitbar(6/16, loadingBar)
fractionRTstructFilename = dir([fractionFolder, 'RTSTRUCT*']).name;
fractionContourData = dicomContours(dicominfo([fractionFolder, fractionRTstructFilename], UseVRHeuristic = false));

%% Extract baseline MR properties

waitbar(7/16, loadingBar, 'Matching baseline data onto fraction image...')

% Number of rows/columns/slices
numberOfBasisRows = basisMRInfo.ImageSize(1);
numberOfBasisColumns = basisMRInfo.ImageSize(2);
numberOfBasisSlices = basisMRInfo.ImageSize(3);

% Voxel dimensions in mm
basisPixelsize = basisMRInfo.PixelSpacings(1);
basisSlicethickness = basisMRInfo.PatientPositions(2, 3) - basisMRInfo.PatientPositions(1, 3);

% Coordinates of center of voxel in first row/column/slice
basisFirstvoxelCenterX = basisMRInfo.PatientPositions(1, 1);
basisFirstvoxelCenterY = basisMRInfo.PatientPositions(1, 2);
basisFirstvoxelCenterZ = basisMRInfo.PatientPositions(1, 3);

% Baseline MR range in mm
basisMRXMin = basisFirstvoxelCenterX - basisPixelsize/2;
basisMRXMax = basisMRXMin + numberOfBasisColumns*basisPixelsize;
basisMRYMin = basisFirstvoxelCenterY - basisPixelsize/2;
basisMRYMax = basisMRYMin + numberOfBasisRows*basisPixelsize;
basisMRZMin = basisFirstvoxelCenterZ - basisSlicethickness/2;
basisMRZMax = basisMRZMin + numberOfBasisSlices*basisSlicethickness;

% Isocenter location in baseline MR
basisIsocROIRownumber = find(string(basisContourData.ROIs{:, 2}) == 'Isocenter');
basisIsocCoordinates = cell2mat(basisContourData.ROIs{basisIsocROIRownumber, 3}{1, 1});
basisIsocVoxel = ...
    [0.5 + (basisIsocCoordinates(2) - basisMRYMin)/basisPixelsize, ...
    0.5 + (basisIsocCoordinates(1) - basisMRXMin)/basisPixelsize, ...
    0.5 + (basisIsocCoordinates(3) - basisMRZMin)/basisSlicethickness];

%% Extract fraction MR properties in similar fashion

numberOfFractionRows = fractionMRInfo.ImageSize(1);
numberOfFractionColumns = fractionMRInfo.ImageSize(2);
numberOfFractionSlices = fractionMRInfo.ImageSize(3);

fractionPixelsize = fractionMRInfo.PixelSpacings(1);
fractionSlicethickness = fractionMRInfo.PatientPositions(2, 3) - fractionMRInfo.PatientPositions(1, 3);
fractionUnitVoxelVolume = fractionPixelsize*fractionPixelsize*fractionSlicethickness;

fractionFirstvoxelCenterX = fractionMRInfo.PatientPositions(1, 1);
fractionFirstvoxelCenterY = fractionMRInfo.PatientPositions(1, 2);
fractionFirstvoxelCenterZ = fractionMRInfo.PatientPositions(1, 3);

fractionMRXMin = fractionFirstvoxelCenterX - fractionPixelsize/2;
fractionMRXMax = fractionMRXMin + numberOfFractionColumns*fractionPixelsize;
fractionMRYMin = fractionFirstvoxelCenterY - fractionPixelsize/2;
fractionMRYMax = fractionMRYMin + numberOfFractionRows*fractionPixelsize;
fractionMRZMin = fractionFirstvoxelCenterZ - fractionSlicethickness/2;
fractionMRZMax = fractionMRZMin + numberOfFractionSlices*fractionSlicethickness;

fractionIsocROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'Isocenter');
fractionIsocCoordinates = cell2mat(fractionContourData.ROIs{fractionIsocROIRownumber, 3}{1, 1});
fractionIsocVoxel = ...
    [0.5 + (fractionIsocCoordinates(2) - fractionMRYMin)/fractionPixelsize, ...
    0.5 + (fractionIsocCoordinates(1) - fractionMRXMin)/fractionPixelsize, ...
    0.5 + (fractionIsocCoordinates(3) - fractionMRZMin)/fractionSlicethickness];

%% Match baseline target delineations and dose distribution onto fraction image

% Shift in x,y,z coordinates (mm) and in number of rows/columns/slices
shiftCoordinates = fractionIsocCoordinates - basisIsocCoordinates;
shiftVoxels = fractionIsocVoxel - basisIsocVoxel;

% Shifting baseline prostate contour point coordinates
basisCTVROIRownumber = find(string(basisContourData.ROIs{:, 2}) == 'CTV');
basisCTVContourpoints = basisContourData.ROIs{basisCTVROIRownumber, 3}{1, 1};
numberOfBasisCTVContours = length(basisCTVContourpoints);
for counter_basiscontours = 1:numberOfBasisCTVContours
    basisCTVContourpointsShifted{counter_basiscontours, 1} = ...
        basisCTVContourpoints{counter_basiscontours, 1} + shiftCoordinates;
end

basisPTVROIRownumber = find(string(basisContourData.ROIs{:, 2}) == 'PTV');
basisPTVContourpoints = basisContourData.ROIs{basisPTVROIRownumber, 3}{1, 1};
numberOfBasisPTVContours = length(basisPTVContourpoints);
for counter_basiscontours = 1:numberOfBasisPTVContours
    basisPTVContourpointsShifted{counter_basiscontours, 1} = ...
        basisPTVContourpoints{counter_basiscontours, 1} + shiftCoordinates;
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

basisDoseDataPercentageShifted = ...
    interp3(matchingFractionColumnNumbers, matchingFractionRowNumbers, matchingFractionSliceNumbers, ...
    basisDoseDataPercentage, basisColumnNumbers, basisRowNumbers, basisSliceNumbers);

%% Extract relevant delineations

waitbar(8/16, loadingBar, 'Performing calculations...')

% Fraction target structures
fractionCTVROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'CTV');
fractionCTVContourpoints = fractionContourData.ROIs{fractionCTVROIRownumber, 3}{1, 1};
fractionPTVROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'PTV');
fractionPTVContourpoints = fractionContourData.ROIs{fractionPTVROIRownumber, 3}{1, 1};

% Fraction OARs in proximity of target structures
fractionBladderROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'BLADDER_IN_RING');
fractionBladderContourpoints = fractionContourData.ROIs{fractionBladderROIRownumber, 3}{1, 1};
fractionRectumROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'RECTUM_IN_RING');
fractionRectumContourpoints = fractionContourData.ROIs{ ...
    fractionRectumROIRownumber, 3}{1, 1};
fractionSmallbowelROIRownumber = find(string(fractionContourData.ROIs{:, 2}) == 'SMALLBOWEL_IN_RI');
fractionSmallbowelContourpoints = fractionContourData.ROIs{fractionSmallbowelROIRownumber, 3}{1, 1};

relevantContourData = deleteContour(basisContourData, 1:basisIsocROIRownumber);
relevantContourData = addContour(relevantContourData, 1, 'CTV_old', basisCTVContourpointsShifted, 'Closed_planar');
relevantContourData = addContour(relevantContourData, 2, 'CTV_new', fractionCTVContourpoints, 'Closed_planar');
relevantContourData = addContour(relevantContourData, 3, 'PTV_old', basisPTVContourpointsShifted, 'Closed_planar');
relevantContourData = addContour(relevantContourData, 4, 'PTV_new', fractionPTVContourpoints, 'Closed_planar');
relevantContourData = addContour(relevantContourData, 5, 'Bladder', fractionBladderContourpoints, 'Closed_planar');
relevantContourData = addContour(relevantContourData, 6, 'Rectum', fractionRectumContourpoints, 'Closed_planar');
% Only include small bowel if it has been delineated
if isempty(fractionSmallbowelContourpoints) == 0
    relevantContourData = addContour(relevantContourData, 7, 'Small bowel', fractionSmallbowelContourpoints, 'Closed_planar');
end

%% Pixelize relevant delineations

fractionMRSpatialInfo = imref3d([numberOfFractionRows, numberOfFractionColumns, numberOfFractionSlices], ...
    [fractionMRXMin + fractionPixelsize/2, fractionMRXMax + fractionPixelsize/2], ...
    [fractionMRYMin + fractionPixelsize/2, fractionMRYMax + fractionPixelsize/2], ...
    [fractionMRZMin + fractionSlicethickness/2, fractionMRZMax + fractionSlicethickness/2]);

oldCTVMask = createMask(relevantContourData, 'CTV_old', fractionMRSpatialInfo);
newCTVMask = createMask(relevantContourData, 'CTV_new', fractionMRSpatialInfo);
oldPTVMask = createMask(relevantContourData, 'PTV_old', fractionMRSpatialInfo);
newPTVMask = createMask(relevantContourData, 'PTV_new', fractionMRSpatialInfo);
bladderMask = createMask(relevantContourData, 'Bladder', fractionMRSpatialInfo);
rectumMask = createMask(relevantContourData, 'Rectum', fractionMRSpatialInfo);
if isempty(fractionSmallbowelContourpoints) == 0
    smallbowelMask = createMask(relevantContourData, 'Small bowel', fractionMRInfo);
% Create empty small bowel mask if it hasn't been delineated
else smallbowelMask = zeros(numberOfFractionRows, numberOfFractionColumns, numberOfFractionSlices);
end

%% Crop fraction MR, shifted dose, and (shifted) delineations

croppingRows = 63:238;
croppingColumns = 32:303;
croppingSlices = 89:200;
numberOfCroppedRows = length(croppingRows);
numberOfCroppedColumns = length(croppingColumns);
numberOfCroppedSlices = length(croppingSlices);

fractionMRDataCropped = fractionMRData(croppingRows, croppingColumns, croppingSlices);
basisDoseDataPercentageShiftedCropped = basisDoseDataPercentageShifted(croppingRows, croppingColumns, croppingSlices);
oldCTVMaskCropped = oldCTVMask(croppingRows, croppingColumns, croppingSlices);
newCTVMaskCropped = newCTVMask(croppingRows, croppingColumns, croppingSlices);
oldPTVMaskCropped = oldPTVMask(croppingRows, croppingColumns, croppingSlices);
newPTVMaskCropped = newPTVMask(croppingRows, croppingColumns, croppingSlices);
bladderMaskCropped = bladderMask(croppingRows, croppingColumns, croppingSlices);
rectumMaskCropped = rectumMask(croppingRows, croppingColumns, croppingSlices);
smallbowelMaskCropped = smallbowelMask(croppingRows, croppingColumns, croppingSlices);

%% Extract target DVH parameters

% First match percent volumes to planned dose in target masks
oldCTVDosePercentages = sort(basisDoseDataPercentageShiftedCropped(oldCTVMaskCropped), 'descend');
numberOfOldCTVVoxels = length(oldCTVDosePercentages);
oldCTVVolume = numberOfOldCTVVoxels*fractionUnitVoxelVolume/1000; %/1000 for conversion from mm^3 to cc
oldCTVAccumulatedPercentageVolume = (1:numberOfOldCTVVoxels)'/numberOfOldCTVVoxels*100;

newCTVDosePercentages = sort(basisDoseDataPercentageShiftedCropped(newCTVMaskCropped), 'descend');
numberOfNewCTVVoxels = length(newCTVDosePercentages);
newCTVVolume = numberOfNewCTVVoxels*fractionUnitVoxelVolume/1000;
newCTVAccumulatedPercentageVolume = (1:numberOfNewCTVVoxels)'/numberOfNewCTVVoxels*100;

oldPTVDosePercentages = sort(basisDoseDataPercentageShiftedCropped(oldPTVMaskCropped), 'descend');
numberOfOldPTVVoxels = length(oldPTVDosePercentages);
oldPTVVolume = numberOfOldPTVVoxels*fractionUnitVoxelVolume/1000;
oldPTVAccumulatedPercentageVolume = (1:numberOfOldPTVVoxels)'/numberOfOldPTVVoxels*100;

newPTVDosePercentages = sort(basisDoseDataPercentageShiftedCropped(newPTVMaskCropped), 'descend');
numberOfNewPTVVoxels = length(newPTVDosePercentages);
newPTVVolume = numberOfNewPTVVoxels*fractionUnitVoxelVolume/1000;
newPTVAccumulatedPercentageVolume = (1:numberOfNewPTVVoxels)'/numberOfNewPTVVoxels*100;

% Then extract DVH parameters
oldCTVD95PlusIndex = find((oldCTVDosePercentages > 95), 1, 'last');
oldCTVD95MinIndex = find((oldCTVDosePercentages < 95), 1);
if isempty(oldCTVD95MinIndex) == 0
    oldCTVV95 = interp1(oldCTVDosePercentages(oldCTVD95PlusIndex:oldCTVD95MinIndex), ...
        oldCTVAccumulatedPercentageVolume(oldCTVD95PlusIndex:oldCTVD95MinIndex), 95);
else oldCTVV95 = 100;
end

oldCTVD107PlusIndex = find((oldCTVDosePercentages > 107), 1, 'last');
oldCTVD107MinIndex = find((oldCTVDosePercentages < 107), 1);
if isempty(oldCTVD107PlusIndex) == 0
    oldCTVV107 = interp1(oldCTVDosePercentages(oldCTVD107PlusIndex:oldCTVD107MinIndex), ...
        oldCTVAccumulatedPercentageVolume(oldCTVD107PlusIndex:oldCTVD107MinIndex), 107);
else oldCTVV107 = 0;
end

newCTVD95PlusIndex = find((newCTVDosePercentages > 95), 1, 'last');
newCTVD95MinIndex = find((newCTVDosePercentages < 95), 1);
if isempty(newCTVD95MinIndex) == 0
    newCTVV95 = interp1(newCTVDosePercentages(newCTVD95PlusIndex:newCTVD95MinIndex), ...
        newCTVAccumulatedPercentageVolume(newCTVD95PlusIndex:newCTVD95MinIndex), 95);
else newCTVV95 = 100;
end

newCTVD107PlusIndex = find((newCTVDosePercentages > 107), 1, 'last');
newCTVD107MinIndex = find((newCTVDosePercentages < 107), 1);
if isempty(newCTVD107PlusIndex) == 0
    newCTVV107 = interp1(newCTVDosePercentages(newCTVD107PlusIndex:newCTVD107MinIndex), ...
        newCTVAccumulatedPercentageVolume(newCTVD107PlusIndex:newCTVD107MinIndex), 107);
else newCTVV107 = 0;
end

% Differences in CTV DVH parameters
differenceCTVVolume = newCTVVolume - oldCTVVolume;
differenceCTVV95 = newCTVV95 - oldCTVV95;
differenceCTVV107 = newCTVV107 - oldCTVV107;


oldPTVD95PlusIndex = find((oldPTVDosePercentages > 95), 1, 'last');
oldPTVD95MinIndex = find((oldPTVDosePercentages < 95), 1);
oldPTVV95 = interp1(oldPTVDosePercentages(oldPTVD95PlusIndex:oldPTVD95MinIndex), ...
    oldPTVAccumulatedPercentageVolume(oldPTVD95PlusIndex:oldPTVD95MinIndex), 95);

oldPTVD107PlusIndex = find((oldPTVDosePercentages > 107), 1, 'last');
oldPTVD107MinIndex = find((oldPTVDosePercentages < 107), 1);
if isempty(oldPTVD107PlusIndex) == 0
    oldPTVV107 = interp1(oldPTVDosePercentages(oldPTVD107PlusIndex:oldPTVD107MinIndex), ...
        oldPTVAccumulatedPercentageVolume(oldPTVD107PlusIndex:oldPTVD107MinIndex), 107);
else oldPTVV107 = 0;
end

newPTVD95PlusIndex = find((newPTVDosePercentages > 95), 1, 'last');
newPTVD95MinIndex = find((newPTVDosePercentages < 95), 1);
newPTVV95 = interp1(newPTVDosePercentages(newPTVD95PlusIndex:newPTVD95MinIndex), ...
    newPTVAccumulatedPercentageVolume(newPTVD95PlusIndex:newPTVD95MinIndex), 95);

newPTVD107PlusIndex = find((newPTVDosePercentages > 107), 1, 'last');
newPTVD107MinIndex = find((newPTVDosePercentages < 107), 1);
if isempty(newPTVD107PlusIndex) == 0
    newPTVV107 = interp1(newPTVDosePercentages(newPTVD107PlusIndex:newPTVD107MinIndex), ...
        newPTVAccumulatedPercentageVolume(newPTVD107PlusIndex:newPTVD107MinIndex), 107);
else newPTVV107 = 0;
end

% Differences in PTV DVH parameters
differencePTVVolume = newPTVVolume - oldPTVVolume;
differencePTVV95 = newPTVV95 - oldPTVV95;
differencePTVV107 = newPTVV107 - oldPTVV107;

%% Set up isodose colourmap

% RGB codes involved colours
maroon = [128, 0, 0]/255;
crimsonred = [220, 20, 60]/255;
neonorange = [255, 95, 31]/255;
orange = [255, 165, 0]/255;
yellow = [255, 255, 50]/255;
green = [0, 255, 0]/255;
mintgreen = [152, 251, 152]/255;
cyan = [0, 255, 255]/255;
navyblue = [0, 0, 128]/255;
magenta = [255, 0, 255]/255;

isodoseLevels = [25; 50; 60; 70; 80; 90; 95; 100; 105; 107];
selectedIsodoseLevels = ones(1, 10);
isodoseColours = [magenta; navyblue; cyan; mintgreen; green; yellow; orange; neonorange; crimsonred; maroon];

% Fill colourmap rows with colours that match the level
isodoseCmap = [];
for counter_isodoses = 1:(length(isodoseLevels) - 1)
    isodoseCmap = [isodoseCmap; ...
        repmat(isodoseColours(counter_isodoses, :), (isodoseLevels(counter_isodoses + 1) - isodoseLevels(counter_isodoses)), 1)];
end
isodoseCmap = [isodoseCmap; repmat(isodoseColours(end, :), (109 - isodoseLevels(end)), 1)];

%% Extract coronal (X, Z), saggital (Y, Z), and transversal (X, Y) planes

corMR = permute(fractionMRDataCropped, [3, 2, 1]);
corDose = permute(basisDoseDataPercentageShiftedCropped, [3, 2, 1]);
corOldCTV = permute(oldCTVMaskCropped, [3, 2, 1]);
corNewCTV = permute(newCTVMaskCropped, [3, 2, 1]);
corOldPTV = permute(oldPTVMaskCropped, [3, 2, 1]);
corNewPTV = permute(newPTVMaskCropped, [3, 2, 1]);
corBladder = permute(bladderMaskCropped, [3, 2, 1]);
corRectum = permute(rectumMaskCropped, [3, 2, 1]);
corSmallbowel = permute(smallbowelMaskCropped, [3, 2, 1]);

sagMR = permute(fractionMRDataCropped, [3, 1, 2]);
sagDose = permute(basisDoseDataPercentageShiftedCropped, [3, 1, 2]);
sagOldCTV = permute(oldCTVMaskCropped, [3, 1, 2]);
sagNewCTV = permute(newCTVMaskCropped, [3, 1, 2]);
sagOldPTV = permute(oldPTVMaskCropped, [3, 1, 2]);
sagNewPTV = permute(newPTVMaskCropped, [3, 1, 2]);
sagBladder = permute(bladderMaskCropped, [3, 1, 2]);
sagRectum = permute(rectumMaskCropped, [3, 1, 2]);
sagSmallbowel = permute(smallbowelMaskCropped, [3, 1, 2]);

transMR = fractionMRDataCropped;
transDose = basisDoseDataPercentageShiftedCropped;
transOldCTV = oldCTVMaskCropped;
transNewCTV = newCTVMaskCropped;
transOldPTV = oldPTVMaskCropped;
transNewPTV = newPTVMaskCropped;
transBladder = bladderMaskCropped;
transRectum = rectumMaskCropped;
transSmallbowel = smallbowelMaskCropped;

%% Compute range of delineations

[~, ~, oldCTVXMin] = ind2sub(size(sagOldCTV), find((sagOldCTV == 1), 1));
[~, ~, oldCTVXMax] = ind2sub(size(sagOldCTV), find((sagOldCTV == 1), 1, 'last'));
[~, ~, oldCTVYMin] = ind2sub(size(corOldCTV), find((corOldCTV == 1), 1));
[~, ~, oldCTVYMax] = ind2sub(size(corOldCTV), find((corOldCTV == 1), 1, 'last'));
[~, ~, oldCTVZMin] = ind2sub(size(transOldCTV), find((transOldCTV == 1), 1));
[~, ~, oldCTVZMax] = ind2sub(size(transOldCTV), find((transOldCTV == 1), 1, 'last'));

[~, ~, newCTVXMin] = ind2sub(size(sagNewCTV), find((sagNewCTV == 1), 1));
[~, ~, newCTVXMax] = ind2sub(size(sagNewCTV), find((sagNewCTV == 1), 1, 'last'));
[~, ~, newCTVYMin] = ind2sub(size(corNewCTV), find((corNewCTV == 1), 1));
[~, ~, newCTVYMax] = ind2sub(size(corNewCTV), find((corNewCTV == 1), 1, 'last'));
[~, ~, newCTVZMin] = ind2sub(size(transNewCTV), find((transNewCTV == 1), 1));
[~, ~, newCTVZMax] = ind2sub(size(transNewCTV), find((transNewCTV == 1), 1, 'last'));

[~, ~, oldPTVXMin] = ind2sub(size(sagOldPTV), find((sagOldPTV == 1), 1));
[~, ~, oldPTVXMax] = ind2sub(size(sagOldPTV), find((sagOldPTV == 1), 1, 'last'));
[~, ~, oldPTVYMin] = ind2sub(size(corOldPTV), find((corOldPTV == 1), 1));
[~, ~, oldPTVYMax] = ind2sub(size(corOldPTV), find((corOldPTV == 1), 1, 'last'));
[~, ~, oldPTVZMin] = ind2sub(size(transOldPTV), find((transOldPTV == 1), 1));
[~, ~, oldPTVZMax] = ind2sub(size(transOldPTV), find((transOldPTV == 1), 1, 'last'));

[~, ~, newPTVXMin] = ind2sub(size(sagNewPTV), find((sagNewPTV == 1), 1));
[~, ~, newPTVXMax] = ind2sub(size(sagNewPTV), find((sagNewPTV == 1), 1, 'last'));
[~, ~, newPTVYMin] = ind2sub(size(corNewPTV), find((corNewPTV == 1), 1));
[~, ~, newPTVYMax] = ind2sub(size(corNewPTV), find((corNewPTV == 1), 1, 'last'));
[~, ~, newPTVZMin] = ind2sub(size(transNewPTV), find((transNewPTV == 1), 1));
[~, ~, newPTVZMax] = ind2sub(size(transNewPTV), find((transNewPTV == 1), 1, 'last'));

[~, ~, bladderXMin] = ind2sub(size(sagBladder), find((sagBladder == 1), 1));
[~, ~, bladderXMax] = ind2sub(size(sagBladder), find((sagBladder == 1), 1, 'last'));
[~, ~, bladderYMin] = ind2sub(size(corBladder), find((corBladder == 1), 1));
[~, ~, bladderYMax] = ind2sub(size(corBladder), find((corBladder == 1), 1, 'last'));
[~, ~, bladderZMin] = ind2sub(size(transBladder), find((transBladder == 1), 1));
[~, ~, bladderZMax] = ind2sub(size(transBladder), find((transBladder == 1), 1, 'last'));

[~, ~, rectumXMin] = ind2sub(size(sagRectum), find((sagRectum == 1), 1));
[~, ~, rectumXMax] = ind2sub(size(sagRectum), find((sagRectum == 1), 1, 'last'));
[~, ~, rectumYMin] = ind2sub(size(corRectum), find((corRectum == 1), 1));
[~, ~, rectumYMax] = ind2sub(size(corRectum), find((corRectum == 1), 1, 'last'));
[~, ~, rectumZMin] = ind2sub(size(transRectum), find((transRectum == 1), 1));
[~, ~, rectumZMax] = ind2sub(size(transRectum), find((transRectum == 1), 1, 'last'));

if isempty(fractionSmallbowelContourpoints) == 0
    [~, ~, smallbowelXMin] = ind2sub(size(sagSmallbowel), find((sagSmallbowel == 1), 1));
    [~, ~, smallbowelXMax] = ind2sub(size(sagSmallbowel), find((sagSmallbowel == 1), 1, 'last'));
    [~, ~, smallbowelYMin] = ind2sub(size(corSmallbowel), find((corSmallbowel == 1), 1));
    [~, ~, smallbowelYMax] = ind2sub(size(corSmallbowel), find((corSmallbowel == 1), 1, 'last'));
    [~, ~, smallbowelZMin] = ind2sub(size(transSmallbowel), find((transSmallbowel == 1), 1));
    [~, ~, smallbowelZMax] = ind2sub(size(transSmallbowel), find((transSmallbowel == 1), 1, 'last'));
% Set small bowel range to zero if it hasn't been delineated 
else [smallbowelXMin, smallbowelXMax, smallbowelYMin, smallbowelYMax, smallbowelZMin, smallbowelZMax] = deal(0);
end

%% Create empty axes in empty uifigure

waitbar(9/16, loadingBar, 'Preparing figures...')

% Set size of coronal/saggital/transversal/DVH plots
corImageratio = numberOfCroppedSlices/numberOfCroppedColumns;
corPlotWidth = 0.44*CZEFullwidth;
corPlotHeight = corPlotWidth*corImageratio;
corPlotLeft = CZEFullwidth/18;
corPlotTop = CZEFullheight - 50;
corPlotBottom = corPlotTop - corPlotHeight;

sagImageratio = numberOfCroppedSlices/numberOfCroppedRows;
sagPlotHeight = corPlotHeight;
sagPlotWidth = sagPlotHeight/sagImageratio;
sagPlotLeft = corPlotLeft + corPlotWidth + 15;
sagPlotBottom = corPlotBottom;

transImageratio = numberOfCroppedRows/numberOfCroppedColumns;
transPlotWidth = corPlotWidth;
transPlotHeight = transPlotWidth*transImageratio;
transPlotLeft = corPlotLeft;
transPlotBottom =  50;

DVHPlotWidth = sagPlotWidth + 10;
DVHPlotHeight = sagPlotHeight;
DVHPlotLeft  = sagPlotLeft;
DVHPlotBottom = sagPlotBottom - DVHPlotHeight - 15;

% Create axes for coronal MR image, dose image, and other plots
axCorMR = uiaxes(f1, 'InnerPosition', [corPlotLeft, corPlotBottom, corPlotWidth, corPlotHeight], ...
    'XLim', [1, numberOfCroppedColumns], 'YLim', [1, numberOfCroppedSlices]);
axCorMR.Visible = 'off';
axCorMR.Toolbar.Visible = 'off';
disableDefaultInteractivity(axCorMR)

axCorDose = uiaxes(f1, 'InnerPosition', [corPlotLeft, corPlotBottom, corPlotWidth, corPlotHeight], ...
    'XLim', [1, numberOfCroppedColumns], 'YLim', [1, numberOfCroppedSlices]);
axCorDose.Visible = 'off';
axCorDose.Toolbar.Visible = 'off';
disableDefaultInteractivity(axCorDose)

axCorOther = uiaxes(f1, 'InnerPosition', [corPlotLeft, corPlotBottom, corPlotWidth, corPlotHeight], ...
    'XLim', [1, numberOfCroppedColumns], 'YLim', [1, numberOfCroppedSlices]);
axCorOther.Visible = 'off';
axCorOther.Toolbar.Visible = 'off';
disableDefaultInteractivity(axCorOther)

linkaxes([axCorMR, axCorDose, axCorOther]);
colormap(axCorMR, 'gray');
clim(axCorDose, [isodoseLevels(1), 109])
colormap(axCorDose, isodoseCmap);

% Saggital axes
axSagMR = uiaxes(f1, 'InnerPosition', [sagPlotLeft, sagPlotBottom, sagPlotWidth, sagPlotHeight], ...
    'XLim', [1, numberOfCroppedRows], 'YLim', [1, numberOfCroppedSlices]);
set(axSagMR, 'XDir', 'reverse');
axSagMR.Visible = 'off';
axSagMR.Toolbar.Visible = 'off';
disableDefaultInteractivity(axSagMR)

axSagDose = uiaxes(f1, 'InnerPosition', [sagPlotLeft, sagPlotBottom, sagPlotWidth, sagPlotHeight], ...
    'XLim', [1, numberOfCroppedRows], 'YLim', [1, numberOfCroppedSlices]);
set(axSagDose, 'XDir', 'reverse');
axSagDose.Visible = 'off';
axSagDose.Toolbar.Visible = 'off';
disableDefaultInteractivity(axSagDose)

axSagOther = uiaxes(f1, 'InnerPosition', [sagPlotLeft, sagPlotBottom, sagPlotWidth, sagPlotHeight], ...
    'XLim', [1, numberOfCroppedRows], 'YLim', [1, numberOfCroppedSlices]);
set(axSagOther, 'XDir', 'reverse');
axSagOther.Visible = 'off';
axSagOther.Toolbar.Visible = 'off';
disableDefaultInteractivity(axSagOther)

linkaxes([axSagMR, axSagDose, axSagOther]);
colormap(axSagMR, 'gray');
clim(axSagDose, [isodoseLevels(1), 109])
colormap(axSagDose, isodoseCmap);

% Transversal axes
axTransMR = uiaxes(f1, 'InnerPosition', [transPlotLeft, transPlotBottom, transPlotWidth, transPlotHeight], ...
    'XLim', [1, numberOfCroppedColumns], 'YLim', [1, numberOfCroppedRows]);
set(axTransMR, 'YDir', 'reverse')
axTransMR.Visible = 'off';
axTransMR.Toolbar.Visible = 'off';
disableDefaultInteractivity(axTransMR)

axTransDose = uiaxes(f1, 'InnerPosition', [transPlotLeft, transPlotBottom, transPlotWidth, transPlotHeight], ...
    'XLim', [1, numberOfCroppedColumns], 'YLim', [1, numberOfCroppedRows]);
set(axTransDose, 'YDir', 'reverse')
axTransDose.Visible = 'off';
axTransDose.Toolbar.Visible = 'off';
disableDefaultInteractivity(axTransDose)

axTransOther = uiaxes(f1, 'InnerPosition', [transPlotLeft, transPlotBottom, transPlotWidth, transPlotHeight], ...
    'XLim', [1, numberOfCroppedColumns], 'YLim', [1, numberOfCroppedRows]);
set(axTransOther, 'YDir', 'reverse');
axTransOther.Visible = 'off';
axTransOther.Toolbar.Visible = 'off';
disableDefaultInteractivity(axTransOther)

linkaxes([axTransMR, axTransDose, axTransOther]);
colormap(axTransMR, 'gray');
clim(axTransDose, [isodoseLevels(1), 109])
colormap(axTransDose, isodoseCmap);

% DVH axis
axDVH = uiaxes(f1, 'OuterPosition', [DVHPlotLeft, DVHPlotBottom, DVHPlotWidth, DVHPlotHeight], ...
    'XLim', [90, 110], 'YLim', [0, 100]);
axDVH.XLabel.Interpreter = 'latex';
axDVH.XLabel.String = 'Dose (\% of prescribed dose)';
axDVH.YLabel.Interpreter = 'latex';
axDVH.YLabel.String = 'Accumulated volume (\%)';
axDVH.Toolbar.Visible = 'off';
disableDefaultInteractivity(axDVH)

%% Create text labels in uifigure

% Title
TitleLabel = uilabel(f1, 'Position', [(CZEFullwidth/2 - 200), (CZEFullheight - 20), 400, 20], ...
    'Text', ['$\textbf{Patient: }$' num2str(selectedPatient) ', $\textbf{fraction: }$ 1'], 'Interpreter', 'latex');
TitleLabel.HorizontalAlignment = 'center';
TitleLabel.VerticalAlignment = 'bottom';

% Coronal plot
XcorLabel = uilabel(f1, 'Position', [(corPlotLeft + corPlotWidth/2 - 10), (corPlotBottom + corPlotHeight + 20), 20, 20], ...
    'Text', 'X', 'Interpreter', 'latex');
XcorLabel.HorizontalAlignment = 'center';
XcorLabel.VerticalAlignment = 'bottom';

LcorLabel = uilabel(f1, 'Position', [(corPlotLeft + corPlotWidth/2 + 125), (corPlotBottom + corPlotHeight), 100, 20], ...
    'Text', 'L $\rightarrow$', 'Interpreter', 'latex');
LcorLabel.HorizontalAlignment = 'left';
LcorLabel.VerticalAlignment = 'bottom';

RcorLabel = uilabel(f1, 'Position', [(corPlotLeft + corPlotWidth/2 - 225), (corPlotBottom + corPlotHeight), 100, 20], ...
    'Text', '$\leftarrow$ R', 'Interpreter', 'latex');
RcorLabel.HorizontalAlignment = 'right';
RcorLabel.VerticalAlignment = 'bottom';

ZcorLabel = uilabel(f1, 'Position', [(corPlotLeft - 45), (corPlotBottom + corPlotHeight/2 - 10), 20, 20], ...
    'Text', 'Z', 'Interpreter', 'latex');
ZcorLabel.HorizontalAlignment = 'right';
ZcorLabel.VerticalAlignment = 'center';

ScorLabel = uilabel(f1, 'Position', [(corPlotLeft - 25), (corPlotBottom + corPlotHeight/2 + 125), 20, 100], ...
    'Text', {'$\uparrow$'; 'S'}, 'Interpreter', 'latex');
ScorLabel.HorizontalAlignment = 'right';
ScorLabel.VerticalAlignment = 'bottom';

IcorLabel = uilabel(f1, 'Position', [(corPlotLeft - 25), (corPlotBottom + corPlotHeight/2 - 225), 20, 100], ...
    'Text', {'I';'$\downarrow$'}, 'Interpreter', 'latex');
IcorLabel.HorizontalAlignment = 'right';
IcorLabel.VerticalAlignment = 'top';

% Saggital plot
YsagLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 10), (sagPlotBottom + sagPlotHeight + 20), 20, 20], ...
    'Text', 'Y', 'Interpreter', 'latex');
YsagLabel.HorizontalAlignment = 'center';
YsagLabel.VerticalAlignment = 'bottom';

AsagLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 + 125), (sagPlotBottom + sagPlotHeight), 100, 20], ...
    'Text', 'A $\rightarrow$', 'Interpreter', 'latex');
AsagLabel.HorizontalAlignment = 'left';
AsagLabel.VerticalAlignment = 'bottom';

PsagLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 225), (sagPlotBottom + sagPlotHeight), 100, 20], ...
    'Text', '$\leftarrow$ P', 'Interpreter', 'latex');
PsagLabel.HorizontalAlignment = 'right';
PsagLabel.VerticalAlignment = 'bottom';

% Transversal plot
YtransLabel = uilabel(f1, 'Position', [(transPlotLeft - 45), (transPlotBottom + transPlotHeight/2 - 10), 20, 20], ...
    'Text', 'Y', 'Interpreter', 'latex');
YtransLabel.HorizontalAlignment = 'right';
YtransLabel.VerticalAlignment = 'center';

AtransLabel = uilabel(f1, 'Position', [(transPlotLeft - 25), (transPlotBottom + transPlotHeight/2 + 125), 20, 100], ...
    'Text', {'$\uparrow$'; 'A'}, 'Interpreter', 'latex');
AtransLabel.HorizontalAlignment = 'right';
AtransLabel.VerticalAlignment = 'bottom';

PtransLabel = uilabel(f1, 'Position', [(transPlotLeft - 25), (transPlotBottom + transPlotHeight/2 - 225), 20, 100], ...
    'Text', {'P'; '$\downarrow$'}, 'Interpreter', 'latex');
PtransLabel.HorizontalAlignment = 'right';
PtransLabel.VerticalAlignment = 'top';

% DVH plot
VolumeLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 110), (DVHPlotBottom + DVHPlotHeight - 115), 160, 100], ...
    'Text', {'$\textbf{Volume:}$'; ...
    ['Old CTV: ' num2str(oldCTVVolume, '%.1f') ' cc']; ...
    ['New CTV: ' num2str(newCTVVolume, '%.1f') ' cc']; ...
    '---------------------------'; ...
    ['Difference: ' num2str(differenceCTVVolume, '%+.1f') ' cc']}, ...
    'Interpreter', 'latex');
VolumeLabel.HorizontalAlignment = 'center';
VolumeLabel.VerticalAlignment = 'top';

V95Label = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 20), (DVHPlotBottom + 40), 160, 200], ...
    'Text', {'$\textbf{V95:}$'; ...
    ['Old CTV: ' num2str(oldCTVV95, '%.1f') '\%']; ...
    ['New CTV: ' num2str(newCTVV95, '%.1f') '\%']; ...
    '-------------------------'; ...
    ['Difference: ' num2str(differenceCTVV95, '%+.1f') '\%']; ''; ''; ...
    ['Old PTV: ' num2str(oldPTVV95, '%.1f') '\%']; ...
    ['New PTV: ' num2str(newPTVV95, '%.1f') '\%']; ...
    '---------------------------'; ...
    ['Difference: ' num2str(differencePTVV95, '%+.1f') '\%']}, ...
    'Interpreter', 'latex');
V95Label.HorizontalAlignment = 'center';
V95Label.VerticalAlignment = 'bottom';

V107Label = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 200), (DVHPlotBottom + 40), 160, 200], ...
    'Text', {'$\textbf{V107:}$'; ...
    ['Old CTV: ' num2str(oldCTVV107, '%.1f') '\%']; ...
    ['New CTV: ' num2str(newCTVV107, '%.1f') '\%']; ...
    '-------------------------'; ...
    ['Difference: ' num2str(differenceCTVV107, '%+.1f') '\%']; ''; ''; ...
    ['Old PTV: ' num2str(oldPTVV107, '%.1f') '\%']; ...
    ['New PTV: ' num2str(newPTVV107, '%.1f') '\%']; ...
    '-------------------------'; ...
    ['Difference: ' num2str(differencePTVV107, '%+.1f') '\%']}, ...
    'Interpreter', 'latex');
V107Label.HorizontalAlignment = 'center';
V107Label.VerticalAlignment = 'bottom';

% Plane selection
PlaneLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 100), (transPlotBottom + 170), 200, 20], ...
    'Text', '$\textbf{Plane selection}$', 'Interpreter', 'latex');
PlaneLabel.HorizontalAlignment = 'center';
PlaneLabel.VerticalAlignment = 'bottom';

XSelectionLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 20), (transPlotBottom + 150), 40, 20], ...
    'Text', 'X:', 'Interpreter', 'latex');
XSelectionLabel.HorizontalAlignment = 'left';
XSelectionLabel.VerticalAlignment = 'bottom';

XSelectionValueLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 20), (transPlotBottom + 150), 40, 20], ...
    'Text', num2str(numberOfCroppedColumns/2), 'Interpreter', 'latex');
XSelectionValueLabel.HorizontalAlignment = 'right';
XSelectionValueLabel.VerticalAlignment = 'bottom';

YSelectionLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 20), (transPlotBottom + 95), 40, 20], ...
    'Text', 'Y:', 'Interpreter', 'latex');
YSelectionLabel.HorizontalAlignment = 'left';
YSelectionLabel.VerticalAlignment = 'bottom';

YSelectionValueLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 20), (transPlotBottom + 95), 40, 20], ...
    'Text', num2str(numberOfCroppedRows/2), 'Interpreter', 'latex');
YSelectionValueLabel.HorizontalAlignment = 'right';
YSelectionValueLabel.VerticalAlignment = 'bottom';

ZSelectionLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 20), (transPlotBottom + 40), 40, 20], ...
    'Text', 'Z:', 'Interpreter', 'latex');
ZSelectionLabel.HorizontalAlignment = 'left';
ZSelectionLabel.VerticalAlignment = 'bottom';

ZSelectionValueLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth/2 - 20), (transPlotBottom + 40), 40, 20], ...
    'Text', num2str(numberOfCroppedSlices/2), 'Interpreter', 'latex');
ZSelectionValueLabel.HorizontalAlignment = 'right';
ZSelectionValueLabel.VerticalAlignment = 'bottom';
 
% Zoom selection
ZoomLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 130), (transPlotBottom + 170), 120, 20], ...
    'Text', '$\textbf{Zoom selection}$', 'Interpreter', 'latex');
ZoomLabel.HorizontalAlignment = 'center';
ZoomLabel.VerticalAlignment = 'bottom';

sagZoomLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 130), (transPlotBottom + 150), 120, 20], ...
    'Text', 'Saggital', 'Interpreter', 'latex');
sagZoomLabel.HorizontalAlignment = 'center';
sagZoomLabel.VerticalAlignment = 'bottom';

corZoomLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 130), (transPlotBottom + 95), 120, 20], ...
    'Text', 'Coronal', 'Interpreter', 'latex');
corZoomLabel.HorizontalAlignment = 'center';
corZoomLabel.VerticalAlignment = 'bottom';

transZoomLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 130), (transPlotBottom + 40), 120, 20], ...
    'Text', 'Transversal', 'Interpreter', 'latex');
transZoomLabel.HorizontalAlignment = 'center';
transZoomLabel.VerticalAlignment = 'bottom';

PlusZoomLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 295), (transPlotBottom + 79), 20, 20], ...
    'Text', '$+$', 'Interpreter', 'latex');
PlusZoomLabel.HorizontalAlignment = 'right';
PlusZoomLabel.VerticalAlignment = 'bottom';

MinZoomLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 60), (transPlotBottom + 77), 20, 20], ...
    'Text', '$-$', 'Interpreter', 'latex');
MinZoomLabel.HorizontalAlignment = 'left';
MinZoomLabel.VerticalAlignment = 'bottom';

% Visibility selection
VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 330), 120, 20], ...
    'Text', '$\textbf{Visibility selection}$', 'Interpreter', 'latex');
VisibilityLabel.HorizontalAlignment = 'center';
VisibilityLabel.VerticalAlignment = 'bottom';

structuresVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 310), 120, 20], ...
    'Text', 'Structures:', 'Interpreter', 'latex');
structuresVisibilityLabel.HorizontalAlignment = 'center';
structuresVisibilityLabel.VerticalAlignment = 'bottom';

oldCTVVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 155), (sagPlotBottom + 290), 60, 20], ...
    'Text', 'Old CTV', 'Interpreter', 'latex');
oldCTVVisibilityLabel.HorizontalAlignment = 'left';
oldCTVVisibilityLabel.VerticalAlignment = 'bottom';

oldPTVVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 155), (sagPlotBottom + 270), 60, 20], ...
    'Text', 'Old PTV', 'Interpreter', 'latex');
oldPTVVisibilityLabel.HorizontalAlignment = 'left';
oldPTVVisibilityLabel.VerticalAlignment = 'bottom';

newCTVVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 270), (sagPlotBottom + 290), 60, 20], ...
    'Text', 'New CTV', 'Interpreter', 'latex');
newCTVVisibilityLabel.HorizontalAlignment = 'left';
newCTVVisibilityLabel.VerticalAlignment = 'bottom';

newPTVVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 270), (sagPlotBottom + 270), 60, 20], ...
    'Text', 'New PTV', 'Interpreter', 'latex');
newPTVVisibilityLabel.HorizontalAlignment = 'left';
newPTVVisibilityLabel.VerticalAlignment = 'bottom';

bladderVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 250), 120, 20], ...
    'Text', 'Bladder', 'Interpreter', 'latex');
bladderVisibilityLabel.HorizontalAlignment = 'center';
bladderVisibilityLabel.VerticalAlignment = 'bottom';

rectumVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 230), 120, 20], ...
    'Text', 'Rectum', 'Interpreter', 'latex');
rectumVisibilityLabel.HorizontalAlignment = 'center';
rectumVisibilityLabel.VerticalAlignment = 'bottom';

smallbowelVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 210), 120, 20], ...
    'Text', 'Small bowel', 'Interpreter', 'latex');
smallbowelVisibilityLabel.HorizontalAlignment = 'center';
smallbowelVisibilityLabel.VerticalAlignment = 'bottom';

isodosesVibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 170), 120, 20], ...
    'Text', 'Isodoses:', 'Interpreter', 'latex');
isodosesVibilityLabel.HorizontalAlignment = 'center';
isodosesVibilityLabel.VerticalAlignment = 'bottom';

isodoseFillingVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 150), 120, 20], ...
    'Text', 'Filling', 'Interpreter', 'latex');
isodoseFillingVisibilityLabel.HorizontalAlignment = 'center';
isodoseFillingVisibilityLabel.VerticalAlignment = 'bottom';

isodose107VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 165), (sagPlotBottom + 130), 50, 20], ...
    'Text', '107\%', 'Interpreter', 'latex');
isodose107VisibilityLabel.HorizontalAlignment = 'right';
isodose107VisibilityLabel.VerticalAlignment = 'bottom';

isodose105VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 165), (sagPlotBottom + 110), 50, 20], ...
    'Text', '105\%', 'Interpreter', 'latex');
isodose105VisibilityLabel.HorizontalAlignment = 'right';
isodose105VisibilityLabel.VerticalAlignment = 'bottom';

isodose100VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 165), (sagPlotBottom + 90), 50, 20], ...
    'Text', '100\%', 'Interpreter', 'latex');
isodose100VisibilityLabel.HorizontalAlignment = 'right';
isodose100VisibilityLabel.VerticalAlignment = 'bottom';

isodose95VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 165), (sagPlotBottom + 70), 50, 20], ...
    'Text', '95\%', 'Interpreter', 'latex');
isodose95VisibilityLabel.HorizontalAlignment = 'right';
isodose95VisibilityLabel.VerticalAlignment = 'bottom';

isodose90VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 165), (sagPlotBottom + 50), 50, 20], ...
    'Text', '90\%', 'Interpreter', 'latex');
isodose90VisibilityLabel.HorizontalAlignment = 'right';
isodose90VisibilityLabel.VerticalAlignment = 'bottom';

isodose80VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 245), (sagPlotBottom + 130), 50, 20], ...
    'Text', '80\%', 'Interpreter', 'latex');
isodose80VisibilityLabel.HorizontalAlignment = 'right';
isodose80VisibilityLabel.VerticalAlignment = 'bottom';

isodose70VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 245), (sagPlotBottom + 110), 50, 20], ...
    'Text', '70\%', 'Interpreter', 'latex');
isodose70VisibilityLabel.HorizontalAlignment = 'right';
isodose70VisibilityLabel.VerticalAlignment = 'bottom';

isodose60VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 245), (sagPlotBottom + 90), 50, 20], ...
    'Text', '60\%', 'Interpreter', 'latex');
isodose60VisibilityLabel.HorizontalAlignment = 'right';
isodose60VisibilityLabel.VerticalAlignment = 'bottom';

isodose50VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 245), (sagPlotBottom + 70), 50, 20], ...
    'Text', '50\%', 'Interpreter', 'latex');
isodose50VisibilityLabel.HorizontalAlignment = 'right';
isodose50VisibilityLabel.VerticalAlignment = 'bottom';

isodose25VisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 245), (sagPlotBottom + 50), 50, 20], ...
    'Text', '25\%', 'Interpreter', 'latex');
isodose25VisibilityLabel.HorizontalAlignment = 'right';
isodose25VisibilityLabel.VerticalAlignment = 'bottom';

XYlinesVisibilityLabel = uilabel(f1, 'Position', [(sagPlotLeft + sagPlotWidth + 170), (sagPlotBottom + 10), 120, 20], ...
    'Text', 'Position indicators', 'Interpreter', 'latex');
XYlinesVisibilityLabel.HorizontalAlignment = 'center';
XYlinesVisibilityLabel.VerticalAlignment = 'bottom';

%% Plot patient data

waitbar(10/16, loadingBar)

% Images have axes starting at left top, plots at left bottom
% Set Y-direction of images to normal to align images with plots (not for transversal plot where Y-direction is reversed)

% Coronal plots
corMRPlot = imagesc(axCorMR, corMR(:, :, numberOfCroppedRows/2));
set(axCorMR, 'YDir', 'normal')

[~, corDosePlot] = contourf(axCorDose, corDose(:, :, numberOfCroppedRows/2), isodoseLevels, ...
    ':', 'FaceAlpha', 0.25, 'LineWidth', 1.5, 'EdgeColor', 'flat');

[~, corOldCTVPlot] = contour(axCorOther, corOldCTV(:, :, numberOfCroppedRows/2), 1, '--', 'Color', 'b', 'LineWidth', 1.5);
hold(axCorOther, 'on')
[~, corOldPTVPlot] = contour(axCorOther, corOldPTV(:, :, numberOfCroppedRows/2), 1, '--', 'Color', 'r', 'LineWidth', 1.5);
[~, corNewCTVPlot] = contour(axCorOther, corNewCTV(:, :, numberOfCroppedRows/2), 1, '-', 'Color', 'b', 'LineWidth', 1.5);
[~, corNewPTVPlot] = contour(axCorOther, corNewPTV(:, :, numberOfCroppedRows/2), 1, '-', 'Color', 'r', 'LineWidth', 1.5);

% Also plot OARs if present in central row
if max(max(corBladder(:, :, numberOfCroppedRows/2))) == 1
    [~, corBladderPlot] = contour(axCorOther, corBladder(:, :, numberOfCroppedRows/2), 1, ...
        '-', 'Color', [254; 242; 0]/255, 'LineWidth', 1.5);
% Create plot, but delete to be plotted data if OAR is not present 
else [~, corBladderPlot] = contour(axCorOther, corBladder(:, :, bladderYMin), 1, ...
        '-', 'Color', [254; 242; 0]/255, 'LineWidth', 1.5);
    corBladderPlot.ZData = [];
end

if max(max(corRectum(:, :, numberOfCroppedRows/2))) == 1
    [~, corRectumPlot] = contour(axCorOther, corRectum(:, :, numberOfCroppedRows/2), 1, ...
        '-', 'Color', [137; 106; 41]/255, 'LineWidth', 1.5);
else [~, corRectumPlot] = contour(axCorOther, corRectum(:, :, rectumYMin), 1, ...
        '-', 'Color', [137; 106; 41]/255, 'LineWidth', 1.5);
    corRectumPlot.ZData = [];
end

if max(max(corSmallbowel(:, :, numberOfCroppedRows/2))) == 1
    [~, corSmallbowelPlot] = contour(axCorOther, corSmallbowel(:, :, numberOfCroppedRows/2), 1, ...
        '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
% Create dummy contour when small bowel hasn't been delineated at all
elseif max(max(max(smallbowelMaskCropped(:, :, :)))) == 0
    [~, corSmallbowelPlot] = contour(axCorOther, eye(5), 1, '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
    corSmallbowelPlot.ZData = [];
else [~, corSmallbowelPlot] = contour(axCorOther, corSmallbowel(:, :, smallbowelYMin), 1, ...
        '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
    corSmallbowelPlot.ZData = [];
end

corXline = xline(axCorOther, numberOfCroppedColumns/2, 'Color', 'y');
corYline = yline(axCorOther, numberOfCroppedSlices/2, 'Color', 'y');
hold(axCorOther, 'off')

waitbar(11/16, loadingBar)

% Saggital plots
sagMRPlot = imagesc(axSagMR, sagMR(:, :, numberOfCroppedColumns/2));
set(axSagMR, 'YDir', 'normal')

[~, sagDosePlot] = contourf(axSagDose, sagDose(:, :, numberOfCroppedColumns/2), isodoseLevels, ...
    ':', 'FaceAlpha', 0.25, 'LineWidth', 1.5, 'EdgeColor', 'flat');

[~, sagOldCTVPlot] = contour(axSagOther, sagOldCTV(:, :, numberOfCroppedColumns/2), 1, '--', 'Color', 'b', 'LineWidth', 1.5);
hold(axSagOther, 'on')
[~, sagOldPTVPlot] = contour(axSagOther, sagOldPTV(:, :, numberOfCroppedColumns/2), 1, '--', 'Color', 'r', 'LineWidth', 1.5);
[~, sagNewCTVPlot] = contour(axSagOther, sagNewCTV(:, :, numberOfCroppedColumns/2), 1, '-', 'Color', 'b', 'LineWidth', 1.5);
[~, sagNewPTVPlot] = contour(axSagOther, sagNewPTV(:, :, numberOfCroppedColumns/2), 1, '-', 'Color', 'r', 'LineWidth', 1.5);

if max(max(sagBladder(:, :, numberOfCroppedColumns/2))) == 1
    [~, sagBladderPlot] = contour(axSagOther, sagBladder(:, :, numberOfCroppedColumns/2), 1, ...
        '-', 'Color', [254; 242; 0]/255, 'LineWidth', 1.5);
else [~, sagBladderPlot] = contour(axSagOther, sagBladder(:, :, bladderXMin), 1, ...
        '-', 'Color', [254; 242; 0]/255, 'LineWidth', 1.5);
    sagBladderPlot.ZData = [];
end

if max(max(sagRectum(:, :, numberOfCroppedColumns/2))) == 1
    [~, sagRectumPlot] = contour(axSagOther, sagRectum(:, :, numberOfCroppedColumns/2), 1, ...
        '-', 'Color', [137; 106; 41]/255, 'LineWidth', 1.5);
else [~, sagRectumPlot] = contour(axSagOther, sagRectum(:, :, rectumXMin), 1, ...
        '-', 'Color', [137; 106; 41]/255, 'LineWidth', 1.5);
    sagRectumPlot.ZData = [];
end

if max(max(sagSmallbowel(:, :, numberOfCroppedColumns/2))) == 1
    [~, sagSmallbowelPlot] = contour(axSagOther, sagSmallbowel(:, :, numberOfCroppedColumns/2), 1, ...
        '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
elseif max(max(max(smallbowelMaskCropped(:, :, :)))) == 0
    [~, sagSmallbowelPlot] = contour(axSagOther, eye(5), 1, '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
    sagSmallbowelPlot.ZData = [];
else [~, sagSmallbowelPlot] = contour(axSagOther, sagSmallbowel(:, :, smallbowelXMin), 1, ...
    '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
end

sagXline = xline(axSagOther, numberOfCroppedRows/2, 'Color', 'y');
sagYline = yline(axSagOther, numberOfCroppedSlices/2, 'Color', 'y');
hold(axSagOther, 'off')

waitbar(12/16, loadingBar)

% Transversal plots
transMRPlot = imagesc(axTransMR, transMR(:, :, numberOfCroppedSlices/2));

[~, transDosePlot] = contourf(axTransDose, transDose(:, :, numberOfCroppedSlices/2), isodoseLevels, ...
    ':', 'FaceAlpha', 0.25, 'LineWidth', 1.5, 'EdgeColor', 'flat');

[~, transOldCTVPlot] = contour(axTransOther, transOldCTV(:, :, numberOfCroppedSlices/2), 1, '--', 'Color', 'b', 'LineWidth', 1.5);
hold(axTransOther, 'on')
[~, transOldPTVPlot] = contour(axTransOther, transOldPTV(:, :, numberOfCroppedSlices/2), 1, '--', 'Color', 'r', 'LineWidth', 1.5);
[~, transNewCTVPlot] = contour(axTransOther, transNewCTV(:, :, numberOfCroppedSlices/2), 1, '-', 'Color', 'b', 'LineWidth', 1.5);
[~, transNewPTVPlot] = contour(axTransOther, transNewPTV(:, :, numberOfCroppedSlices/2), 1, '-', 'Color', 'r', 'LineWidth', 1.5);

if max(max(transBladder(:, :, numberOfCroppedSlices/2))) == 1
    [~, transBladderPlot] = contour(axTransOther, transBladder(:, :, numberOfCroppedSlices/2), 1, ...
        '-', 'Color', [254; 242; 0]/255, 'LineWidth', 1.5);
else [~, transBladderPlot] = contour(axTransOther, transBladder(:, :, bladderZMin), 1, ...
        '-', 'Color', [254; 242; 0]/255, 'LineWidth', 1.5);
    transBladderPlot.ZData = [];
end

if max(max(transRectum(:, :, numberOfCroppedSlices/2))) == 1
    [~, transRectumPlot] = contour(axTransOther, transRectum(:, :, numberOfCroppedSlices/2), 1, ...
        '-', 'Color', [137; 106; 41]/255, 'LineWidth', 1.5);
else [~, transRectumPlot] = contour(axTransOther, transRectum(:, :, rectumZMin), 1, ...
        '-', 'Color', [137; 106; 41]/255, 'LineWidth', 1.5);
    transRectumPlot.ZData = [];
end

if max(max(transSmallbowel(:, :, numberOfCroppedSlices/2))) == 1
    [~, transSmallbowelPlot] = contour(axTransOther, transSmallbowel(:, :, numberOfCroppedSlices/2), 1, ...
        '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
elseif max(max(max(smallbowelMaskCropped(:, :, :)))) == 0
    [~, transSmallbowelPlot] = contour(axTransOther, eye(5), 1, '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
    transSmallbowelPlot.ZData = [];
else [~, transSmallbowelPlot] = contour(axTransOther, transSmallbowel(:, :, smallbowelZMin), 1, ...
        '-', 'Color', [188; 102; 0]/255, 'LineWidth', 1.5);
    transSmallbowelPlot.ZData = [];
end

transXline = xline(axTransOther, numberOfCroppedColumns/2, 'Color', 'y');
transYline = yline(axTransOther, numberOfCroppedRows/2, 'Color', 'y');
hold(axTransOther, 'off')

waitbar(13/16, loadingBar)

% Common isodose colorbar
cbarDose = colorbar(axCorDose, 'Units', 'pixels', 'Position', ...
    [(sagPlotLeft + sagPlotWidth + 15), sagPlotBottom, 50, sagPlotHeight]);
cbarDose.Ticks = isodoseLevels;
cbarDose.Label.Interpreter = 'latex';
cbarDose.Label.String = 'Dose (\% of prescribed dose)';
cbarDose.FontSize = 10;

% DVH curve
plot(axDVH, oldCTVDosePercentages, oldCTVAccumulatedPercentageVolume, ...
    '--', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'Old CTV')
hold(axDVH, 'on')
plot(axDVH, oldPTVDosePercentages, oldPTVAccumulatedPercentageVolume, ...
    '--', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'Old PTV')
plot(axDVH, newCTVDosePercentages, newCTVAccumulatedPercentageVolume, ...
    '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'New CTV')
plot(axDVH, newPTVDosePercentages, newPTVAccumulatedPercentageVolume, ...
    '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'New PTV')
xline(axDVH, 95, '--', 'HandleVisibility', 'off')
xline(axDVH, 107, '--', 'HandleVisibility', 'off')
hold(axDVH, 'off')
axDVH.FontSize = 13;
legend(axDVH, 'Location', 'NorthEast', 'FontSize', 10);
   
%% Create sliders, buttons, and checkboxes

waitbar(14/16, loadingBar, 'Finishing...')

% X selection slider and buttons
XSelection = uislider(f1, 'Position', [sagPlotLeft + 30, ...
    (transPlotBottom + 140), (sagPlotWidth - 60), 3], ...
    'Limits', [1, numberOfCroppedColumns], 'Value', numberOfCroppedColumns/2, ...
    'MajorTicks', [1, numberOfCroppedColumns/2, numberOfCroppedColumns], ...
    'MajorTickLabels', {'1', '', num2str(numberOfCroppedColumns)}, ...
    'MinorTicks', 4:4:numberOfCroppedColumns, ...
    'ValueChangingFcn', @(src, evt) changingXSelection(src, evt, ...
    XSelectionValueLabel, sagMR, sagMRPlot, sagDose, sagDosePlot, ...
    sagOldCTV, sagOldCTVPlot, oldCTVXMin, oldCTVXMax, ...
    sagOldPTV, sagOldPTVPlot, oldPTVXMin, oldPTVXMax, ...
    sagNewCTV, sagNewCTVPlot, newCTVXMin, newCTVXMax, ...
    sagNewPTV, sagNewPTVPlot, newPTVXMin, newPTVXMax, ...
    sagBladder, sagBladderPlot, bladderXMin, bladderXMax, ...
    sagRectum, sagRectumPlot, rectumXMin, rectumXMax, ...
    sagSmallbowel, sagSmallbowelPlot, smallbowelXMin, smallbowelXMax, ...
    corXline, transXline));

XDown = uibutton(f1, 'Position', [sagPlotLeft, ...
    (transPlotBottom + 133), 20, 20], 'Text', '<', ...
    'ButtonPushedFcn', @(src, evt) pushXDown(src, evt, ...
    XSelection, XSelectionValueLabel, sagMR, sagMRPlot, sagDose, sagDosePlot, ...
    sagOldCTV, sagOldCTVPlot, oldCTVXMin, oldCTVXMax, ...
    sagOldPTV, sagOldPTVPlot, oldPTVXMin, oldPTVXMax, ...
    sagNewCTV, sagNewCTVPlot, newCTVXMin, newCTVXMax, ...
    sagNewPTV, sagNewPTVPlot, newPTVXMin, newPTVXMax, ...
    sagBladder, sagBladderPlot, bladderXMin, bladderXMax, ...
    sagRectum, sagRectumPlot, rectumXMin, rectumXMax, ...
    sagSmallbowel, sagSmallbowelPlot, smallbowelXMin, smallbowelXMax, ...
    corXline, transXline));

XUp = uibutton(f1, 'Position', [(sagPlotLeft + sagPlotWidth - 20), ...
    (transPlotBottom + 133), 20, 20], 'Text', '>', ...
    'ButtonPushedFcn', @(src, evt) pushXUp(src, evt, ...
    numberOfCroppedColumns, XSelection, XSelectionValueLabel, sagMR, sagMRPlot, ...
    sagDose, sagDosePlot, sagOldCTV, sagOldCTVPlot, oldCTVXMin, oldCTVXMax, ...
    sagOldPTV, sagOldPTVPlot, oldPTVXMin, oldPTVXMax, ...
    sagNewCTV, sagNewCTVPlot, newCTVXMin, newCTVXMax, ...
    sagNewPTV, sagNewPTVPlot, newPTVXMin, newPTVXMax, ...
    sagBladder, sagBladderPlot, bladderXMin, bladderXMax, ...
    sagRectum, sagRectumPlot, rectumXMin, rectumXMax, ...
    sagSmallbowel, sagSmallbowelPlot, smallbowelXMin, smallbowelXMax, ...
    corXline, transXline));

% Y selection slider and buttons
YSelection = uislider(f1, 'Position', [(sagPlotLeft + 30), ...
    (transPlotBottom + 85), (sagPlotWidth - 60), 3], ...
    'Limits', [1, numberOfCroppedRows], 'Value', numberOfCroppedRows/2, ...
    'MajorTicks', [1, numberOfCroppedRows/2, numberOfCroppedRows], ...
    'MajorTickLabels', {'1', '', num2str(numberOfCroppedRows)}, ...
    'MinorTicks', 4:4:numberOfCroppedRows, ...
    'ValueChangingFcn', @(src, evt) changingYSelection(src, evt, ...
    YSelectionValueLabel, corMR, corMRPlot, corDose, corDosePlot, ...
    corOldCTV, corOldCTVPlot, oldCTVYMin, oldCTVYMax, ...
    corOldPTV, corOldPTVPlot, oldPTVYMin, oldPTVYMax, ...
    corNewCTV, corNewCTVPlot, newCTVYMin, newCTVYMax, ...
    corNewPTV, corNewPTVPlot, newPTVYMin, newPTVYMax, ...
    corBladder, corBladderPlot, bladderYMin, bladderYMax, ...
    corRectum, corRectumPlot, rectumYMin, rectumYMax, ...
    corSmallbowel, corSmallbowelPlot, smallbowelYMin, smallbowelYMax, ...
    sagXline, transYline));

YDown = uibutton(f1, 'Position', [sagPlotLeft, ...
    (transPlotBottom + 78), 20, 20], 'Text', '<', ...
    'ButtonPushedFcn', @(src, evt) pushYDown(src, evt, ...
    YSelection, YSelectionValueLabel, corMR, corMRPlot, corDose, corDosePlot, ...
    corOldCTV, corOldCTVPlot, oldCTVYMin, oldCTVYMax, ...
    corOldPTV, corOldPTVPlot, oldPTVYMin, oldPTVYMax, ...
    corNewCTV, corNewCTVPlot, newCTVYMin, newCTVYMax, ...
    corNewPTV, corNewPTVPlot, newPTVYMin, newPTVYMax, ...
    corBladder, corBladderPlot, bladderYMin, bladderYMax, ...
    corRectum, corRectumPlot, rectumYMin, rectumYMax, ...
    corSmallbowel, corSmallbowelPlot, smallbowelYMin, smallbowelYMax, ...
    sagXline, transYline));

YUp = uibutton(f1, 'Position', [(sagPlotLeft + sagPlotWidth - 20), ...
    (transPlotBottom + 78), 20, 20], 'Text', '>', ...
    'ButtonPushedFcn', @(src, evt) pushYUp(src, evt, ...
    numberOfCroppedRows, YSelection, YSelectionValueLabel, corMR, corMRPlot, ...
    corDose, corDosePlot, corOldCTV, corOldCTVPlot, oldCTVYMin, oldCTVYMax, ...
    corOldPTV, corOldPTVPlot, oldPTVYMin, oldPTVYMax, ...
    corNewCTV, corNewCTVPlot, newCTVYMin, newCTVYMax, ...
    corNewPTV, corNewPTVPlot, newPTVYMin, newPTVYMax, ...
    corBladder, corBladderPlot, bladderYMin, bladderYMax, ...
    corRectum, corRectumPlot, rectumYMin, rectumYMax, ...
    corSmallbowel, corSmallbowelPlot, smallbowelYMin, smallbowelYMax, ...
    sagXline, transYline));

% Z selection slider and buttons
ZSelection = uislider(f1, 'Position', [(sagPlotLeft + 30), ...
    (transPlotBottom + 30), (sagPlotWidth - 60), 3], ...
    'Limits', [1, numberOfCroppedSlices], 'Value', numberOfCroppedSlices/2, ...
    'MajorTicks', [1, numberOfCroppedSlices/2, numberOfCroppedSlices], ...
    'MajorTickLabels', {'1', '', num2str(numberOfCroppedSlices)}, ...
    'MinorTicks', 4:4:numberOfCroppedSlices, ...
    'ValueChangingFcn', @(src, evt) changingZSelection(src, evt, ...
    ZSelectionValueLabel, transMR, transMRPlot, transDose, transDosePlot, ...
    transOldCTV, transOldCTVPlot, oldCTVZMin, oldCTVZMax, ...
    transOldPTV, transOldPTVPlot, oldPTVZMin, oldPTVZMax, ...
    transNewCTV, transNewCTVPlot, newCTVZMin, newCTVZMax, ...
    transNewPTV, transNewPTVPlot, newPTVZMin, newPTVZMax, ...
    transBladder, transBladderPlot, bladderZMin, bladderZMax, ...
    transRectum, transRectumPlot, rectumZMin, rectumZMax, ...
    transSmallbowel, transSmallbowelPlot, smallbowelZMin, smallbowelZMax, ...
    corYline, sagYline));

ZDown = uibutton(f1, 'Position', [sagPlotLeft, ...
    (transPlotBottom + 23), 20, 20], 'Text', '<', ...
    'ButtonPushedFcn', @(src, evt) pushZDown(src, evt, ...
    ZSelection, ZSelectionValueLabel, transMR, transMRPlot, transDose, transDosePlot, ...
    transOldCTV, transOldCTVPlot, oldCTVZMin, oldCTVZMax, ...
    transOldPTV, transOldPTVPlot, oldPTVZMin, oldPTVZMax, ...
    transNewCTV, transNewCTVPlot, newCTVZMin, newCTVZMax, ...
    transNewPTV, transNewPTVPlot, newPTVZMin, newPTVZMax, ...
    transBladder, transBladderPlot, bladderZMin, bladderZMax, ...
    transRectum, transRectumPlot, rectumZMin, rectumZMax, ...
    transSmallbowel, transSmallbowelPlot, smallbowelZMin, smallbowelZMax, ...
    corYline, sagYline));

ZUp = uibutton(f1, 'Position', [(sagPlotLeft + sagPlotWidth - 20), ...
    (transPlotBottom + 23), 20, 20], 'Text', '>', ...
    'ButtonPushedFcn', @(src, evt) pushZUp(src, evt, ...
    numberOfCroppedSlices, ZSelection, ZSelectionValueLabel, transMR, transMRPlot, ...
    transDose, transDosePlot, transOldCTV, transOldCTVPlot, oldCTVZMin, oldCTVZMax, ...
    transOldPTV, transOldPTVPlot, oldPTVZMin, oldPTVZMax, ...
    transNewCTV, transNewCTVPlot, newCTVZMin, newCTVZMax, ...
    transNewPTV, transNewPTVPlot, newPTVZMin, newPTVZMax, ...
    transBladder, transBladderPlot, bladderZMin, bladderZMax, ...
    transRectum, transRectumPlot, rectumZMin, rectumZMax, ...
    transSmallbowel, transSmallbowelPlot, smallbowelZMin, smallbowelZMax, ...
    corYline, sagYline));

% Zoom sliders
waitbar(15/16, loadingBar)

sagZoom = uislider(f1, 'Position', [sagPlotLeft + sagPlotWidth + 90, ...
    (transPlotBottom + 140), 200, 3], ...
    'Limits', [-0.2, 0.2], 'Value', 0, 'MajorTicks', [], 'MinorTicks', [], ...
    'ValueChangingFcn', @(src, evt) changingSagZoom(src, evt, axSagMR), ...
    'ValueChangedFCn', @(src, evt) changedSagZoom(src, evt));

corZoom = uislider(f1, 'Position', [sagPlotLeft + sagPlotWidth + 90, ...
    (transPlotBottom + 85), 200, 3], ...
    'Limits', [-0.2, 0.2], 'Value', 0, 'MajorTicks', [], 'MinorTicks', [], ...
    'ValueChangingFcn', @(src, evt) changingCorZoom(src, evt, axCorMR), ...
    'ValueChangedFCn', @(src, evt) changedCorZoom(src, evt));

transZoom = uislider(f1, 'Position', [sagPlotLeft + sagPlotWidth + 90, ...
    (transPlotBottom + 30), 200, 3], ...
    'Limits', [-0.2, 0.2], 'Value', 0, 'MajorTicks', [], 'MinorTicks', [], ...
    'ValueChangingFcn', @(src, evt) changingTransZoom(src, evt, axTransMR), ...
    'ValueChangedFCn', @(src, evt) changedTransZoom(src, evt));

% Visibility checkboxes
waitbar(1, loadingBar)

oldCTVVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 130, ...
    (sagPlotBottom + 290 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedOldCTVVisibility(src, evt, ...
    corOldCTVPlot, sagOldCTVPlot, transOldCTVPlot));

oldPTVVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 130, ...
    (sagPlotBottom + 270 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedOldPTVVisibility(src, evt, ...
    corOldPTVPlot, sagOldPTVPlot, transOldPTVPlot));

newCTVVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 290 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedNewCTVVisibility(src, evt, ...
    corNewCTVPlot, sagNewCTVPlot, transNewCTVPlot));

newPTVVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 270 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedNewPTVVisibility(src, evt, ...
    corNewPTVPlot, sagNewPTVPlot, transNewPTVPlot));

bladderVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 185, ...
    (sagPlotBottom + 250 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedBladderVisibility(src, evt, ...
    corBladderPlot, sagBladderPlot, transBladderPlot));

rectumVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 185, ...
    (sagPlotBottom + 230 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedRectumVisibility(src, evt, ...
    corRectumPlot, sagRectumPlot, transRectumPlot));

smallbowelVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 170, ...
    (sagPlotBottom + 210 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedSmallbowelVisibility(src, evt, ...
    corSmallbowelPlot, sagSmallbowelPlot, transSmallbowelPlot));

isodoseFillingVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 185, ...
    (sagPlotBottom + 150 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedIsodoseFillingVisibility(src, evt, ...
    corDosePlot, sagDosePlot, transDosePlot));

isodose107Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 165, ...
    (sagPlotBottom + 130 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose105Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 165, ...
    (sagPlotBottom + 110 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose100Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 165, ...
    (sagPlotBottom + 90 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose95Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 165, ...
    (sagPlotBottom + 70 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose90Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 165, ...
    (sagPlotBottom + 50 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose80Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 130 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose70Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 110 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose60Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 90 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose50Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 70 - 1), 20, 20], 'Text', '', 'Value', 1);

isodose25Visibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 245, ...
    (sagPlotBottom + 50 - 1), 20, 20], 'Text', '', 'Value', 1);

% Set one common callback function for all isodose visibility checkboxes
[isodose107Visibility.ValueChangedFcn, isodose105Visibility.ValueChangedFcn, ...
    isodose100Visibility.ValueChangedFcn, isodose95Visibility.ValueChangedFcn, ...
    isodose90Visibility.ValueChangedFcn, isodose80Visibility.ValueChangedFcn, ...
    isodose70Visibility.ValueChangedFcn, isodose60Visibility.ValueChangedFcn, ...
    isodose50Visibility.ValueChangedFcn, isodose25Visibility.ValueChangedFcn] ...
    = deal(@(src, evt) changedIsodoseVisibility(src, evt, corDosePlot, sagDosePlot, transDosePlot, ...
    isodoseLevels, isodose107Visibility, isodose105Visibility, isodose100Visibility, ...
    isodose95Visibility, isodose90Visibility, isodose80Visibility, isodose70Visibility, ...
    isodose60Visibility, isodose50Visibility, isodose25Visibility));

XYlinesVisibility = uicheckbox(f1, 'Position', [sagPlotLeft + sagPlotWidth + 150, ...
    (sagPlotBottom + 10 - 1), 20, 20], 'Text', '', 'Value', 1, ...
    'ValueChangedFcn', @(src, evt) changedXYlinesVisibility(src, evt, ...
    corXline, corYline, sagXline, sagYline, transXline, transYline));

delete(loadingBar)

%% Functions to update plots with sliders, buttons, and checkboxes

% X selection functions
function changingXSelection(src, evt, ...
    XSelectionValueLabel, sagMR, sagMRPlot, sagDose, sagDosePlot, ...
    sagOldCTV, sagOldCTVPlot, oldCTVXMin, oldCTVXMax, ...
    sagOldPTV, sagOldPTVPlot, oldPTVXMin, oldPTVXMax, ...
    sagNewCTV, sagNewCTVPlot, newCTVXMin, newCTVXMax, ...
    sagNewPTV, sagNewPTVPlot, newPTVXMin, newPTVXMax, ...
    sagBladder, sagBladderPlot, bladderXMin, bladderXMax, ...
    sagRectum, sagRectumPlot, rectumXMin, rectumXMax, ...
    sagSmallbowel, sagSmallbowelPlot, smallbowelXMin, smallbowelXMax, ...
    corXline, transXline)

    X = round(evt.Value);
    
    % Update slider value label
    XSelectionValueLabel.Text = num2str(X);
         
    % Show new MR and dose images
    sagMRPlot.CData = sagMR(:, :, X);
    sagDosePlot.ZData = sagDose(:, :, X);
    
    % Show new delineations if applicable
    if X >= oldCTVXMin && X <= oldCTVXMax
        sagOldCTVPlot.ZData = double(sagOldCTV(:, :, X));
    else sagOldCTVPlot.ZData = [];
    end
    
    if X >= oldPTVXMin && X <= oldPTVXMax
        sagOldPTVPlot.ZData = double(sagOldPTV(:, :, X));
    else sagOldPTVPlot.ZData = [];
    end

    if X >= newCTVXMin && X <= newCTVXMax
        sagNewCTVPlot.ZData = double(sagNewCTV(:, :, X));
    else sagNewCTVPlot.ZData = [];
    end
    
    if X >= newPTVXMin && X <= newPTVXMax
        sagNewPTVPlot.ZData = double(sagNewPTV(:, :, X));
    else sagNewPTVPlot.ZData = [];
    end

    if X >= bladderXMin && X <= bladderXMax
            sagBladderPlot.ZData = double(sagBladder(:, :, X));
    else sagBladderPlot.ZData = [];
    end
    
    if X >= rectumXMin && X <= rectumXMax
            sagRectumPlot.ZData = double(sagRectum(:, :, X));
    else sagRectumPlot.ZData = [];
    end
    
    if X >= smallbowelXMin && X <= smallbowelXMax
            sagSmallbowelPlot.ZData = double(sagSmallbowel(:, :, X));
    else sagSmallbowelPlot.ZData = [];
    end

    % Update x/y-lines
    corXline.Value = X;
    transXline.Value = X;
end

function pushXDown(src, evt, XSelection, ...
    XSelectionValueLabel, sagMR, sagMRPlot, sagDose, sagDosePlot, ...
    sagOldCTV, sagOldCTVPlot, oldCTVXMin, oldCTVXMax, ...
    sagOldPTV, sagOldPTVPlot, oldPTVXMin, oldPTVXMax, ...
    sagNewCTV, sagNewCTVPlot, newCTVXMin, newCTVXMax, ...
    sagNewPTV, sagNewPTVPlot, newPTVXMin, newPTVXMax, ...
    sagBladder, sagBladderPlot, bladderXMin, bladderXMax, ...
    sagRectum, sagRectumPlot, rectumXMin, rectumXMax, ...
    sagSmallbowel, sagSmallbowelPlot, smallbowelXMin, smallbowelXMax, ...
    corXline, transXline)
    
    % Check current slider value
    oldX = round(XSelection.Value);
    
    % Update everything as long as lower boundary hasn't been reached
    if oldX > 1
        newX = oldX - 1;
        
        XSelection.Value = newX;
        XSelectionValueLabel.Text = num2str(newX);
        sagMRPlot.CData = sagMR(:, :, newX);
        sagDosePlot.ZData = sagDose(:, :, newX);
    
        if newX >= oldCTVXMin && newX <= oldCTVXMax
            sagOldCTVPlot.ZData = double(sagOldCTV(:, :, newX));
        else sagOldCTVPlot.ZData = [];
        end
    
        if newX >= oldPTVXMin && newX <= oldPTVXMax
            sagOldPTVPlot.ZData = double(sagOldPTV(:, :, newX));
        else sagOldPTVPlot.ZData = [];
        end

        if newX >= newCTVXMin && newX <= newCTVXMax
            sagNewCTVPlot.ZData = double(sagNewCTV(:, :, newX));
        else sagNewCTVPlot.ZData = [];
        end
    
        if newX >= newPTVXMin && newX <= newPTVXMax
            sagNewPTVPlot.ZData = double(sagNewPTV(:, :, newX));
        else sagNewPTVPlot.ZData = [];
        end
        
        if newX >= bladderXMin && newX <= bladderXMax
            sagBladderPlot.ZData = double(sagBladder(:, :, newX));
        else sagBladderPlot.ZData = [];
        end
    
        if newX >= rectumXMin && newX <= rectumXMax
            sagRectumPlot.ZData = double(sagRectum(:, :, newX));
        else sagRectumPlot.ZData = [];
        end
        
        if newX >= smallbowelXMin && newX <= smallbowelXMax
            sagSmallbowelPlot.ZData = double(sagSmallbowel(:, :, newX));
        else sagSmallbowelPlot.ZData = [];
        end

        corXline.Value = newX;
        transXline.Value = newX;
    end
end

function pushXUp(src, evt, numberOfCroppedColumns, XSelection, ...
    XSelectionValueLabel, sagMR, sagMRPlot, sagDose, sagDosePlot, ...
    sagOldCTV, sagOldCTVPlot, oldCTVXMin, oldCTVXMax, ...
    sagOldPTV, sagOldPTVPlot, oldPTVXMin, oldPTVXMax, ...
    sagNewCTV, sagNewCTVPlot, newCTVXMin, newCTVXMax, ...
    sagNewPTV, sagNewPTVPlot, newPTVXMin, newPTVXMax, ...
    sagBladder, sagBladderPlot, bladderXMin, bladderXMax, ...
    sagRectum, sagRectumPlot, rectumXMin, rectumXMax, ...
    sagSmallbowel, sagSmallbowelPlot, smallbowelXMin, smallbowelXMax, ...
    corXline, transXline)
    
    % Check current slider value
    oldX = round(XSelection.Value);
    
    % Update everything as long as upper boundary hasn't been reached
    if oldX < numberOfCroppedColumns
        newX = oldX + 1;
        
        XSelection.Value = newX;
        XSelectionValueLabel.Text = num2str(newX);
        sagMRPlot.CData = sagMR(:, :, newX);
        sagDosePlot.ZData = sagDose(:, :, newX);
    
        if newX >= oldCTVXMin && newX <= oldCTVXMax
            sagOldCTVPlot.ZData = double(sagOldCTV(:, :, newX));
        else sagOldCTVPlot.ZData = [];
        end
    
        if newX >= oldPTVXMin && newX <= oldPTVXMax
            sagOldPTVPlot.ZData = double(sagOldPTV(:, :, newX));
        else sagOldPTVPlot.ZData = [];
        end
        
        if newX >= newCTVXMin && newX <= newCTVXMax
            sagNewCTVPlot.ZData = double(sagNewCTV(:, :, newX));
        else sagNewCTVPlot.ZData = [];
        end
    
        if newX >= newPTVXMin && newX <= newPTVXMax
            sagNewPTVPlot.ZData = double(sagNewPTV(:, :, newX));
        else sagNewPTVPlot.ZData = [];
        end

        if newX >= bladderXMin && newX <= bladderXMax
            sagBladderPlot.ZData = double(sagBladder(:, :, newX));
        else sagBladderPlot.ZData = [];
        end
    
        if newX >= rectumXMin && newX <= rectumXMax
            sagRectumPlot.ZData = double(sagRectum(:, :, newX));
        else sagRectumPlot.ZData = [];
        end
        
        if newX >= smallbowelXMin && newX <= smallbowelXMax
            sagSmallbowelPlot.ZData = double(sagSmallbowel(:, :, newX));
        else sagSmallbowelPlot.ZData = [];
        end

        corXline.Value = newX;
        transXline.Value = newX;
    end
end

% Y selection functions
function changingYSelection(src, evt, ...
    YSelectionValueLabel, corMR, corMRPlot, corDose, corDosePlot, ...
    corOldCTV, corOldCTVPlot, oldCTVYMin, oldCTVYMax, ...
    corOldPTV, corOldPTVPlot, oldPTVYMin, oldPTVYMax, ...
    corNewCTV, corNewCTVPlot, newCTVYMin, newCTVYMax, ...
    corNewPTV, corNewPTVPlot, newPTVYMin, newPTVYMax, ...
    corBladder, corBladderPlot, bladderYMin, bladderYMax, ...
    corRectum, corRectumPlot, rectumYMin, rectumYMax, ...
    corSmallbowel, corSmallbowelPlot, smallbowelYMin, smallbowelYMax, ...
    sagXline, transYline)

    Y = round(evt.Value);
    
    % Update slider value label
    YSelectionValueLabel.Text = num2str(Y);
       
    % Show new MR and dose image
    corMRPlot.CData = corMR(:, :, Y);
    corDosePlot.ZData = corDose(:, :, Y);

    % Show new delineations if applicable
    if Y >= oldCTVYMin && Y <= oldCTVYMax
        corOldCTVPlot.ZData = double(corOldCTV(:, :, Y));
    else corOldCTVPlot.ZData = [];
    end
    
    if Y >= oldPTVYMin && Y <= oldPTVYMax
        corOldPTVPlot.ZData = double(corOldPTV(:, :, Y));
    else corOldPTVPlot.ZData = [];
    end
    
    if Y >= newCTVYMin && Y <= newCTVYMax
        corNewCTVPlot.ZData = double(corNewCTV(:, :, Y));
    else corNewCTVPlot.ZData = [];
    end
    
    if Y >= newPTVYMin && Y <= newPTVYMax
        corNewPTVPlot.ZData = double(corNewPTV(:, :, Y));
    else corNewPTVPlot.ZData = [];
    end
    
    if Y >= bladderYMin && Y <= bladderYMax
            corBladderPlot.ZData = double(corBladder(:, :, Y));
    else corBladderPlot.ZData = [];
    end
    
    if Y >= rectumYMin && Y <= rectumYMax
            corRectumPlot.ZData = double(corRectum(:, :, Y));
    else corRectumPlot.ZData = [];
    end
    
    if Y >= smallbowelYMin && Y <= smallbowelYMax
            corSmallbowelPlot.ZData = double(corSmallbowel(:, :, Y));
    else corSmallbowelPlot.ZData = [];
    end

    % Update x/y-lines
    sagXline.Value = Y;
    transYline.Value = Y;
end

function pushYDown(src, evt, YSelection, ...
    YSelectionValueLabel, corMR, corMRPlot, corDose, corDosePlot, ...
    corOldCTV, corOldCTVPlot, oldCTVYMin, oldCTVYMax, ...
    corOldPTV, corOldPTVPlot, oldPTVYMin, oldPTVYMax, ...
    corNewCTV, corNewCTVPlot, newCTVYMin, newCTVYMax, ...
    corNewPTV, corNewPTVPlot, newPTVYMin, newPTVYMax, ...
    corBladder, corBladderPlot, bladderYMin, bladderYMax, ...
    corRectum, corRectumPlot, rectumYMin, rectumYMax, ...
    corSmallbowel, corSmallbowelPlot, smallbowelYMin, smallbowelYMax, ...
    sagXline, transYline)
    
    % Check current slider value
    oldY = round(YSelection.Value);
    
    % Update everything as long as lower boundary hasn't been reached
    if oldY > 1
        newY = oldY - 1;
        
        YSelection.Value = newY;
        YSelectionValueLabel.Text = num2str(newY);
        corMRPlot.CData = corMR(:, :, newY);
        corDosePlot.ZData = corDose(:, :, newY);

        if newY >= oldCTVYMin && newY <= oldCTVYMax
            corOldCTVPlot.ZData = double(corOldCTV(:, :, newY));
        else corOldCTVPlot.ZData = [];
        end
    
        if newY >= oldPTVYMin && newY <= oldPTVYMax
            corOldPTVPlot.ZData = double(corOldPTV(:, :, newY));
        else corOldPTVPlot.ZData = [];
        end
        
        if newY >= newCTVYMin && newY <= newCTVYMax
            corNewCTVPlot.ZData = double(corNewCTV(:, :, newY));
        else corNewCTVPlot.ZData = [];
        end
    
        if newY >= newPTVYMin && newY <= newPTVYMax
            corNewPTVPlot.ZData = double(corNewPTV(:, :, newY));
        else corNewPTVPlot.ZData = [];
        end

        if newY >= bladderYMin && newY <= bladderYMax
            corBladderPlot.ZData = double(corBladder(:, :, newY));
        else corBladderPlot.ZData = [];
        end
    
        if newY >= rectumYMin && newY <= rectumYMax
            corRectumPlot.ZData = double(corRectum(:, :, newY));
        else corRectumPlot.ZData = [];
        end
        
        if newY >= smallbowelYMin && newY <= smallbowelYMax
            corSmallbowelPlot.ZData = double(corSmallbowel(:, :, newY));
        else corSmallbowelPlot.ZData = [];
        end

        sagXline.Value = newY;
        transYline.Value = newY;
    end
end

function pushYUp(src, evt, numberOfCroppedRows, YSelection, ...
    YSelectionValueLabel, corMR, corMRPlot, corDose, corDosePlot, ...
    corOldCTV, corOldCTVPlot, oldCTVYMin, oldCTVYMax, ...
    corOldPTV, corOldPTVPlot, oldPTVYMin, oldPTVYMax, ...
    corNewCTV, corNewCTVPlot, newCTVYMin, newCTVYMax, ...
    corNewPTV, corNewPTVPlot, newPTVYMin, newPTVYMax, ...
    corBladder, corBladderPlot, bladderYMin, bladderYMax, ...
    corRectum, corRectumPlot, rectumYMin, rectumYMax, ...
    corSmallbowel, corSmallbowelPlot, smallbowelYMin, smallbowelYMax, ...
    sagXline, transYline)
    
    % Check current slider value
    oldY = round(YSelection.Value);
    
    % Update everything as long as upper boundary hasn't been reached
    if oldY < numberOfCroppedRows
        newY = oldY + 1;
        
        YSelection.Value = newY;
        YSelectionValueLabel.Text = num2str(newY);
        corMRPlot.CData = corMR(:, :, newY);
        corDosePlot.ZData = corDose(:, :, newY);

        if newY >= oldCTVYMin && newY <= oldCTVYMax
            corOldCTVPlot.ZData = double(corOldCTV(:, :, newY));
        else corOldCTVPlot.ZData = [];
        end
    
        if newY >= oldPTVYMin && newY <= oldPTVYMax
            corOldPTVPlot.ZData = double(corOldPTV(:, :, newY));
        else corOldPTVPlot.ZData = [];
        end

        if newY >= newCTVYMin && newY <= newCTVYMax
            corNewCTVPlot.ZData = double(corNewCTV(:, :, newY));
        else corNewCTVPlot.ZData = [];
        end
    
        if newY >= newPTVYMin && newY <= newPTVYMax
            corNewPTVPlot.ZData = double(corNewPTV(:, :, newY));
        else corNewPTVPlot.ZData = [];
        end

        if newY >= bladderYMin && newY <= bladderYMax
            corBladderPlot.ZData = double(corBladder(:, :, newY));
        else corBladderPlot.ZData = [];
        end
    
        if newY >= rectumYMin && newY <= rectumYMax
            corRectumPlot.ZData = double(corRectum(:, :, newY));
        else corRectumPlot.ZData = [];
        end

        if newY >= smallbowelYMin && newY <= smallbowelYMax
            corSmallbowelPlot.ZData = double(corSmallbowel(:, :, newY));
        else corSmallbowelPlot.ZData = [];
        end
            
        sagXline.Value = newY;
        transYline.Value = newY;
    end
end

% Z selection functions
function changingZSelection(src, evt, ...
    ZSelectionValueLabel, transMR, transMRPlot, transDose, transDosePlot, ...
    transOldCTV, transOldCTVPlot, oldCTVZMin, oldCTVZMax, ...
    transOldPTV, transOldPTVPlot, oldPTVZMin, oldPTVZMax, ...
    transNewCTV, transNewCTVPlot, newCTVZMin, newCTVZMax, ...
    transNewPTV, transNewPTVPlot, newPTVZMin, newPTVZMax, ...
    transBladder, transBladderPlot, bladderZMin, bladderZMax, ...
    transRectum, transRectumPlot, rectumZMin, rectumZMax, ...
    transSmallbowel, transSmallbowelPlot, smallbowelZMin, smallbowelZMax, ...
    corYline, sagYline)

    Z = round(evt.Value);
    
    % Update slider value label
    ZSelectionValueLabel.Text = num2str(Z);
    
    % Show new MR and dose image
    transMRPlot.CData = transMR(:, :, Z);
    transDosePlot.ZData = transDose(:, :, Z);
    
    % Show new delineations if applicable
    if Z >= oldCTVZMin && Z <= oldCTVZMax
        transOldCTVPlot.ZData = double(transOldCTV(:, :, Z));
    else transOldCTVPlot.ZData = [];
    end
    
    if Z >= oldPTVZMin && Z <= oldPTVZMax
        transOldPTVPlot.ZData = double(transOldPTV(:, :, Z));
    else transOldPTVPlot.ZData = [];
    end

    if Z >= newCTVZMin && Z <= newCTVZMax
        transNewCTVPlot.ZData = double(transNewCTV(:, :, Z));
    else transNewCTVPlot.ZData = [];
    end
    
    if Z >= newPTVZMin && Z <= newPTVZMax
        transNewPTVPlot.ZData = double(transNewPTV(:, :, Z));
    else transNewPTVPlot.ZData = [];
    end

    if Z >= bladderZMin && Z <= bladderZMax
            transBladderPlot.ZData = double(transBladder(:, :, Z));
    else transBladderPlot.ZData = [];
    end
        
    if Z >= rectumZMin && Z <= rectumZMax
            transRectumPlot.ZData = double(transRectum(:, :, Z));
    else transRectumPlot.ZData = [];
    end
    
    if Z >= smallbowelZMin && Z <= smallbowelZMax
            transSmallbowelPlot.ZData = double(transSmallbowel(:, :, Z));
    else transSmallbowelPlot.ZData = [];
    end

    % Update x/y-lines
    corYline.Value = Z;
    sagYline.Value = Z;
end

function pushZDown(src, evt, ZSelection, ...
    ZSelectionValueLabel, transMR, transMRPlot, transDose, transDosePlot, ...
    transOldCTV, transOldCTVPlot, oldCTVZMin, oldCTVZMax, ...
    transOldPTV, transOldPTVPlot, oldPTVZMin, oldPTVZMax, ...
    transNewCTV, transNewCTVPlot, newCTVZMin, newCTVZMax, ...
    transNewPTV, transNewPTVPlot, newPTVZMin, newPTVZMax, ...
    transBladder, transBladderPlot, bladderZMin, bladderZMax, ...
    transRectum, transRectumPlot, rectumZMin, rectumZMax, ...
    transSmallbowel, transSmallbowelPlot, smallbowelZMin, smallbowelZMax, ...
    corYline, sagYline)
    
    % Check current slider value
    oldZ = round(ZSelection.Value);
    
    % Update everything as long as lower boundary hasn't been reached
    if oldZ > 1
        newZ = oldZ - 1;
        
        ZSelection.Value = newZ;
        ZSelectionValueLabel.Text = num2str(newZ);
        transMRPlot.CData = transMR(:, :, newZ);
        transDosePlot.ZData = transDose(:, :, newZ);
    
        if newZ >= oldCTVZMin && newZ <= oldCTVZMax
            transOldCTVPlot.ZData = double(transOldCTV(:, :, newZ));
        else transOldCTVPlot.ZData = [];
        end
    
        if newZ >= oldPTVZMin && newZ <= oldPTVZMax
            transOldPTVPlot.ZData = double(transOldPTV(:, :, newZ));
        else transOldPTVPlot.ZData = [];
        end

        if newZ >= newCTVZMin && newZ <= newCTVZMax
            transNewCTVPlot.ZData = double(transNewCTV(:, :, newZ));
        else transNewCTVPlot.ZData = [];
        end
    
        if newZ >= newPTVZMin && newZ <= newPTVZMax
            transNewPTVPlot.ZData = double(transNewPTV(:, :, newZ));
        else transNewPTVPlot.ZData = [];
        end

        if newZ >= bladderZMin && newZ <= bladderZMax
            transBladderPlot.ZData = double(transBladder(:, :, newZ));
        else transBladderPlot.ZData = [];
        end
        
        if newZ >= rectumZMin && newZ <= rectumZMax
            transRectumPlot.ZData = double(transRectum(:, :, newZ));
        else transRectumPlot.ZData = [];
        end
        
        if newZ >= smallbowelZMin && newZ <= smallbowelZMax
            transSmallbowelPlot.ZData = double(transSmallbowel(:, :, newZ));
        else transSmallbowelPlot.ZData = [];
        end

        corYline.Value = newZ;
        sagYline.Value = newZ;
    end
end

function pushZUp(src, evt, numberOfCroppedSlices, ZSelection, ...
    ZSelectionValueLabel, transMR, transMRPlot, transDose, transDosePlot, ...
    transOldCTV, transOldCTVPlot, oldCTVZMin, oldCTVZMax, ...
    transOldPTV, transOldPTVPlot, oldPTVZMin, oldPTVZMax, ...
    transNewCTV, transNewCTVPlot, newCTVZMin, newCTVZMax, ...
    transNewPTV, transNewPTVPlot, newPTVZMin, newPTVZMax, ...
    transBladder, transBladderPlot, bladderZMin, bladderZMax, ...
    transRectum, transRectumPlot, rectumZMin, rectumZMax, ...
    transSmallbowel, transSmallbowelPlot, smallbowelZMin, smallbowelZMax, ...
    corYline, sagYline)
    
    % Check current slider value
    oldZ = round(ZSelection.Value);
    
    % Update everything as long as upper boundary hasn't been reached
    if oldZ < numberOfCroppedSlices
        newZ = oldZ + 1;
        
        ZSelection.Value = newZ;      
        ZSelectionValueLabel.Text = num2str(newZ);
        
        ZSelection.Value = newZ;
        ZSelectionValueLabel.Text = num2str(newZ);
        transMRPlot.CData = transMR(:, :, newZ);
        transDosePlot.ZData = transDose(:, :, newZ);
    
        if newZ >= oldCTVZMin && newZ <= oldCTVZMax
            transOldCTVPlot.ZData = double(transOldCTV(:, :, newZ));
        else transOldCTVPlot.ZData = [];
        end
    
        if newZ >= oldPTVZMin && newZ <= oldPTVZMax
            transOldPTVPlot.ZData = double(transOldPTV(:, :, newZ));
        else transOldPTVPlot.ZData = [];
        end

        if newZ >= newCTVZMin && newZ <= newCTVZMax
            transNewCTVPlot.ZData = double(transNewCTV(:, :, newZ));
        else transNewCTVPlot.ZData = [];
        end
    
        if newZ >= newPTVZMin && newZ <= newPTVZMax
            transNewPTVPlot.ZData = double(transNewPTV(:, :, newZ));
        else transNewPTVPlot.ZData = [];
        end
        
        if newZ >= bladderZMin && newZ <= bladderZMax
            transBladderPlot.ZData = double(transBladder(:, :, newZ));
        else transBladderPlot.ZData = [];
        end
        
        if newZ >= rectumZMin && newZ <= rectumZMax
            transRectumPlot.ZData = double(transRectum(:, :, newZ));
        else transRectumPlot.ZData = [];
        end

        if newZ >= smallbowelZMin && newZ <= smallbowelZMax
            transSmallbowelPlot.ZData = double(transSmallbowel(:, :, newZ));
        else transSmallbowelPlot.ZData = [];
        end

        corYline.Value = newZ;
        sagYline.Value = newZ;
    end
end

% Zoom functions
function changingSagZoom(src, evt, axSagMR)
    zoomFactor = 10^(evt.Value);
    zoom(axSagMR, zoomFactor)
end

function changedSagZoom(src, evt)
    src.Value = 0;
end

function changingCorZoom(src, evt, axCorMR)
    zoomFactor = 10^(evt.Value);
    zoom(axCorMR, zoomFactor)
end

function changedCorZoom(src, evt)
    src.Value = 0;
end

function changingTransZoom(src, evt, axTransMR)
    zoomFactor = 10^(evt.Value);
    zoom(axTransMR, zoomFactor)
end

function changedTransZoom(src, evt)
    src.Value = 0;
end

% Visibility functions
function changedOldCTVVisibility(src, evt, corOldCTVPlot, sagOldCTVPlot, transOldCTVPlot)
    if evt.Value == 1
        corOldCTVPlot.Visible = 'on';
        sagOldCTVPlot.Visible = 'on';
        transOldCTVPlot.Visible = 'on';
    else
        corOldCTVPlot.Visible = 'off';
        sagOldCTVPlot.Visible = 'off';
        transOldCTVPlot.Visible = 'off';
    end
end

function changedOldPTVVisibility(src, evt, corOldPTVPlot, sagOldPTVPlot, transOldPTVPlot)
    if evt.Value == 1
        corOldPTVPlot.Visible = 'on';
        sagOldPTVPlot.Visible = 'on';
        transOldPTVPlot.Visible = 'on';
    else
        corOldPTVPlot.Visible = 'off';
        sagOldPTVPlot.Visible = 'off';
        transOldPTVPlot.Visible = 'off';
    end
end

function changedNewCTVVisibility(src, evt, corNewCTVPlot, sagNewCTVPlot, transNewCTVPlot)
    if evt.Value == 1
        corNewCTVPlot.Visible = 'on';
        sagNewCTVPlot.Visible = 'on';
        transNewCTVPlot.Visible = 'on';
    else
        corNewCTVPlot.Visible = 'off';
        sagNewCTVPlot.Visible = 'off';
        transNewCTVPlot.Visible = 'off';
    end
end

function changedNewPTVVisibility(src, evt, corNewPTVPlot, sagNewPTVPlot, transNewPTVPlot)
    if evt.Value == 1
        corNewPTVPlot.Visible = 'on';
        sagNewPTVPlot.Visible = 'on';
        transNewPTVPlot.Visible = 'on';
    else
        corNewPTVPlot.Visible = 'off';
        sagNewPTVPlot.Visible = 'off';
        transNewPTVPlot.Visible = 'off';
    end
end

function changedBladderVisibility(src, evt, corBladderPlot, sagBladderPlot, transBladderPlot)
    if evt.Value == 1
        corBladderPlot.Visible = 'on';
        sagBladderPlot.Visible = 'on';
        transBladderPlot.Visible = 'on';
    else
        corBladderPlot.Visible = 'off';
        sagBladderPlot.Visible = 'off';
        transBladderPlot.Visible = 'off';
    end
end

function changedRectumVisibility(src, evt, corRectumPlot, sagRectumPlot, transRectumPlot)
    if evt.Value == 1
        corRectumPlot.Visible = 'on';
        sagRectumPlot.Visible = 'on';
        transRectumPlot.Visible = 'on';
    else
        corRectumPlot.Visible = 'off';
        sagRectumPlot.Visible = 'off';
        transRectumPlot.Visible = 'off';
    end
end

function changedSmallbowelVisibility(src, evt, corSmallbowelPlot, sagSmallbowelPlot, transSmallbowelPlot)
    if evt.Value == 1
        corSmallbowelPlot.Visible = 'on';
        sagSmallbowelPlot.Visible = 'on';
        transSmallbowelPlot.Visible = 'on';
    else
        corSmallbowelPlot.Visible = 'off';
        sagSmallbowelPlot.Visible = 'off';
        transSmallbowelPlot.Visible = 'off';
    end
end

function changedIsodoseFillingVisibility(src, evt, corDosePlot, sagDosePlot, transDosePlot)
    if evt.Value == 1
        corDosePlot.FaceAlpha = 0.25;
        sagDosePlot.FaceAlpha = 0.25;
        transDosePlot.FaceAlpha = 0.25;
    else
        corDosePlot.FaceAlpha = 0;
        sagDosePlot.FaceAlpha = 0;
        transDosePlot.FaceAlpha = 0;
    end
end

function changedIsodoseVisibility(src, evt, corDosePlot, sagDosePlot, transDosePlot, ...
    isodoseLevels, isodose107Visibility, isodose105Visibility, isodose100Visibility, ...
    isodose95Visibility, isodose90Visibility, isodose80Visibility, isodose70Visibility, ...
    isodose60Visibility, isodose50Visibility, isodose25Visibility)
    if isodose107Visibility.Value == 1
        isodoseLevels(10) = 107;
    else
        isodoseLevels(10) = 1000;
    end
    
    if isodose105Visibility.Value == 1
        isodoseLevels(9) = 105;
    else
        isodoseLevels(9) = 1000;
    end
    
    if isodose100Visibility.Value == 1
        isodoseLevels(8) = 100;
    else
        isodoseLevels(8) = 1000;
    end
    
    if isodose95Visibility.Value == 1
        isodoseLevels(7) = 95;
    else
        isodoseLevels(7) = 1000;
    end

    if isodose90Visibility.Value == 1
        isodoseLevels(6) = 90;
    else
        isodoseLevels(6) = 1000;
    end

    if isodose80Visibility.Value == 1
        isodoseLevels(5) = 80;
    else
        isodoseLevels(5) = 1000;
    end

    if isodose70Visibility.Value == 1
        isodoseLevels(4) = 70;
    else
        isodoseLevels(4) = 1000;
    end

    if isodose60Visibility.Value == 1
        isodoseLevels(3) = 60;
    else
        isodoseLevels(3) = 1000;
    end

    if isodose50Visibility.Value == 1
        isodoseLevels(2) = 50;
    else
        isodoseLevels(2) = 1000;
    end

    if isodose25Visibility.Value == 1
        isodoseLevels(1) = 25;
    else
        isodoseLevels(1) = 1000;
    end

    corDosePlot.LevelList = isodoseLevels;
    sagDosePlot.LevelList = isodoseLevels;
    transDosePlot.LevelList = isodoseLevels;
end

function changedXYlinesVisibility(src, evt, corXline, corYline, ...
    sagXline, sagYline, transXline, transYline)
    if evt.Value == 1
        corXline.Visible = 'on';
        corYline.Visible = 'on';
        sagXline.Visible = 'on';
        sagYline.Visible = 'on';
        transXline.Visible = 'on';
        transYline.Visible = 'on';
    else
        corXline.Visible = 'off';
        corYline.Visible = 'off';
        sagXline.Visible = 'off';
        sagYline.Visible = 'off';
        transXline.Visible = 'off';
        transYline.Visible = 'off';
    end
end