clear all;
close all;

img=imread("Dataset\8.png");

denoisedImg = imgaussfilt3(img);

% HSV range for each colors
Black.hueMin = 0.000;
Black.hueMax = 1.000;
Black.satMin = 0.000;
Black.satMax = 0.200;
Black.valMin = 0.000;
Black.valMax = 0.500;

White.hueMin = 0.000;
White.hueMax = 1.000;
White.satMin = 0.000;
White.satMax = 0.200;
White.valMin = 0.700;
White.valMax = 1.000;

Blue.hueMin = 0.500;
Blue.hueMax = 0.700;
Blue.satMin = 0.500;
Blue.satMax = 1.000;
Blue.valMin = 0.500;
Blue.valMax = 1.000;

Red.hueMin = 0.865;
Red.hueMax = 0.075;
Red.satMin = 0.200;
Red.satMax = 1.000;
Red.valMin = 0.300;
Red.valMax = 1.000;

Yellow.hueMin = 0.100;
Yellow.hueMax = 0.200;
Yellow.satMin = 0.550;
Yellow.satMax = 1.000;
Yellow.valMin = 0.550;
Yellow.valMax = 1.000;

Green.hueMin = 0.275;
Green.hueMax = 0.550;
Green.satMin = 0.500;
Green.satMax = 1.000;
Green.valMin = 0.250;
Green.valMax = 1.000;

% Extract based on color
[blackBW, blackMaskRGBImg] = createMask(denoisedImg, Black);
[whiteBW, whiteMaskRGBImg] = createMask(denoisedImg, White);
[blueBW, blueMaskRGBImg] = createMask(denoisedImg, Blue);
[redBW, redMaskRGBImg] = createMask(denoisedImg, Red);
[yellowBW, yellowMaskRGBImg] = createMask(denoisedImg, Yellow);
[greenBW, greenMaskRGBImg] = createMask(denoisedImg, Green);

rowCount = 4;
figure;
subplot(rowCount,5,1);
imshow(img);

% subplot(rowCount,6,1);
% imshow(blackMaskRGBImg);title('black');
% subplot(rowCount,6,7);
% imshow(blackBW);
% 
% subplot(rowCount,6,2);
% imshow(whiteMaskRGBImg);title('white');
% subplot(rowCount,6,8);
% imshow(whiteBW);
% 
% subplot(rowCount,6,3);
% imshow(blueMaskRGBImg);title('blue');
% subplot(rowCount,6,9);
% imshow(blueBW);
% 
% subplot(rowCount,6,4);
% imshow(redMaskRGBImg);title('red');
% subplot(rowCount,6,10);
% imshow(redBW);
% 
% subplot(rowCount,6,5);
% imshow(yellowMaskRGBImg);title('yellow');
% subplot(rowCount,6,11);
% imshow(yellowBW);
% 
% subplot(rowCount,6,6);
% imshow(greenMaskRGBImg);title('green');
% subplot(rowCount,6,12);
% imshow(greenBW);

%
blackBW = morphBW(blackBW);
subplot(rowCount,5,2);
imshow(blackBW);title('black');

blackFilterBW = findCCandFilter(blackBW);
subplot(rowCount,5,3);
imshow(blackFilterBW);title('black filtered');

whiteBW = morphBW(whiteBW);
subplot(rowCount,5,4);
imshow(whiteBW);title('white');

whiteFilterBW = findCCandFilter(whiteBW);
subplot(rowCount,5,5);
imshow(whiteFilterBW);title('white filtered');


blueBW = morphRGBY(blueBW);
subplot(rowCount,5,6);
imshow(blueBW);title('blue');

blueFilterBW = findCCandFilter(blueBW);
subplot(rowCount,5,7);
imshow(blueFilterBW);title('blue filtered');


redBW = morphRGBY(redBW);
subplot(rowCount,5,8);
imshow(redBW);title('red');

redFilterBW = findCCandFilter(redBW);
subplot(rowCount,5,9);
imshow(redFilterBW);title('red filtered');


yellowBW = morphRGBY(yellowBW);
subplot(rowCount,5,10);
imshow(yellowBW);title('yellow');

yellowFilterBW = findCCandFilter(yellowBW);
subplot(rowCount,5,11);
imshow(yellowFilterBW);title('yellow filtered');


greenBW = morphRGBY(greenBW);
subplot(rowCount,5,12);
imshow(greenBW);title('green');

greenFilterBW = findCCandFilter(greenBW);
subplot(rowCount,5,13);
imshow(greenFilterBW);title('green filtered');


% combine all colored CCs
combinedcc = combineBRGY(blueFilterBW, redFilterBW, yellowFilterBW, greenFilterBW);

whitecc = bwconncomp(whiteFilterBW);
blackcc = bwconncomp(blackFilterBW);

% combine white and black CCs
whiteBlackCC = combineCC(whitecc, blackcc);

%%
% edge detection then morphology
gray = rgb2gray(img);

% subplot(rowCount, 5, 1);
% imshow(gray);title('from rgb');

nlmGray = imnlmfilt(gray);
% subplot(rowCount, 5, 2);
% imshow(nlmGray);title('nlmGray');


edgeDetect= edge(nlmGray, 'sobel', 'vertical');
% subplot(rowCount, 5, 3 );
% imshow(edgeDetect);title('edgeDetect');

dilateEdge = imdilate(edgeDetect, strel('rectangle', [5 10]));
% subplot(rowCount, 5, 4 );
% imshow(dilateEdge);title('dilateEdge');

closeEdge = imclose(dilateEdge, strel('rectangle', [5 15]));
% subplot(rowCount, 5, 5 );
% imshow(closeEdge);title('closeEdge');

openEdge = imopen(closeEdge, strel('rectangle', [17 30]));
% subplot(rowCount, 5, 9 );
% imshow(openEdge);title('openEdge');

dilateEdge2 = imdilate(openEdge, strel('rectangle', [1 10]));
% subplot(rowCount, 5, 10 );
% imshow(dilateEdge2);title('dilateEdge2');

closeEdge2 = imclose(dilateEdge2, strel('rectangle', [2 15]));
subplot(rowCount, 5, 16 );
imshow(closeEdge2);title('closeEdge2');


final = findCCandFilter(closeEdge2);
subplot(rowCount, 5, 17);
imshow(final);title('edge final');

edgeCC = bwconncomp(final);


%%
plateFound = false;
maxRatio = 0;

fprintf("checking color");
overlapRatio = checkOverlap(combinedcc, edgeCC)
if ~isempty(overlapRatio)
    [maxRatioVal, maxRatioIdx] = max(overlapRatio, [], 'all');
    [combinedccIdx, edgeCCIdx] = ind2sub(size(overlapRatio), maxRatioIdx);

    if maxRatioVal > 0.6 % overlap threshold
        plateFound = true;
        plateCC = combinedcc;
        plateIdx = combinedccIdx;
    elseif maxRatioVal > 0
        maxRatio = maxRatioVal;
        plateCC = combinedcc;
        plateIdx = combinedccIdx;
    end
else % check with black and white
    fprintf("checking color with BW");
    overlapRatio = checkOverlap(combinedcc, whiteBlackCC)
    if ~isempty(overlapRatio)
    [maxRatioVal, maxRatioIdx] = max(overlapRatio, [], 'all');
    [combinedccIdx, whiteBlackCC] = ind2sub(size(overlapRatio), maxRatioIdx);

        if maxRatioVal > 0.6 % overlap threshold
            plateFound = true;
            plateCC = combinedcc;
            plateIdx = combinedccIdx;
        elseif maxRatioVal > 0
            maxRatio = maxRatioVal;
            plateCC = combinedcc;
            plateIdx = combinedccIdx;
        end
    end
end

% coloured not found proceed to check edge with white
if ~plateFound
    fprintf("checking white");
    overlapRatio = checkOverlap(whitecc, edgeCC)
    if ~isempty(overlapRatio)
        [maxRatioVal, maxRatioIdx] = max(overlapRatio, [], 'all');
        [whiteccIdx, edgeCCIdx] = ind2sub(size(overlapRatio), maxRatioIdx);

        if maxRatioVal > 0.5 % overlap threshold
            plateFound = true;
            plateCC = whitecc; % segment based on white cc regions
            plateIdx = whiteccIdx; % segment based on white cc regions
        elseif maxRatioVal > maxRatio
            maxRatio = maxRatioVal;
            plateCC = whitecc; % segment based on white cc regions
            plateIdx = whiteccIdx; % segment based on white cc regions
        end
    end
end

% white not found proceed to check edge with black
if ~plateFound
    fprintf("checking black");
    overlapRatio = checkOverlap(blackcc, edgeCC)
    if ~isempty(overlapRatio)
        [maxRatioVal, maxRatioIdx] = max(overlapRatio, [], 'all');
        [blackccIdx, edgeCCIdx] = ind2sub(size(overlapRatio), maxRatioIdx);
        
        if maxRatioVal > 0.4 % overlap threshold
            plateFound = true;
            plateCC = edgeCC; % segment based on edge cc regions
            plateIdx = edgeCCIdx;
        elseif maxRatioVal > maxRatio
            maxRatio = maxRatioVal;
            plateCC = edgeCC; % segment based on edge cc regions
            plateIdx = edgeCCIdx;
        end
    end
end

% none found proceed to edge 
if ~plateFound
    if edgeCC.NumObjects == 1 % (gets if there is only 1 cc)
        plateFound = true;
        plateCC = edgeCC;
        plateIdx = 1;
    end
end

% if still no found and got max overlap cc
if plateFound == false && maxRatio > 0
    plateFound = true;
    fprintf("Found based on small overlap");
end


if plateFound
     plateProp = regionprops(plateCC);
     plate = imcrop(img, plateProp(plateIdx).BoundingBox);
     figure;
     imshow(plate); title('plate');
     %imwrite(plate,'plate.png');
     
    
     plateSize = size(plate);
     roiBB = [1 1 plateSize(2), plateSize(1)];
     roi = roiBB;

      % change to true to manually draw own roi, 
      % false to automatically select the whole image as roi
     if (false)
         figure
         imshow(plate);
         title("Draw region of interest (ROI)");

         % Draw a region of interest
          h = imrect
    
          % Evaluate OCR within ROI
          roi = h.getPosition;
     end
      
      ocrI = evaluateOCRTraining(plate, roi);

      % Show results
      figure
      imshow(ocrI);
      title("Recognised characters");
else
    fprintf("Unable to detect license plate");
end




% extract based on HSV threshold
function [BW,maskedRGBImage] = createMask(RGB, dataStruct)

% Convert RGB image to chosen color space
I = rgb2hsv(RGB);

% Define thresholds for hue
hueMin = dataStruct.hueMin;
hueMax = dataStruct.hueMax;

% Define thresholds for saturation
satMin = dataStruct.satMin;
satMax = dataStruct.satMax;

% Define thresholds for value
valMin = dataStruct.valMin;
valMax = dataStruct.valMax;

% Create mask based on chosen thresholds

if hueMin > hueMax
    sliderBW = ( (I(:,:,1) >= hueMin) | (I(:,:,1) <= hueMax) ) & ...
    (I(:,:,2) >= satMin ) & (I(:,:,2) <= satMax) & ...
    (I(:,:,3) >= valMin ) & (I(:,:,3) <= valMax);
else
    sliderBW = ( (I(:,:,1) >= hueMin) & (I(:,:,1) <= hueMax) ) & ...
    (I(:,:,2) >= satMin ) & (I(:,:,2) <= satMax) & ...
    (I(:,:,3) >= valMin ) & (I(:,:,3) <= valMax);
end

BW = bwareaopen(sliderBW, 25); %remove blobs smaller than 25 px

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end

function morphedBW = morphRGBY(BW) % morphology operation for colored image
    morphedBW = imdilate(BW, strel('line', 5, 0));
    morphedBW = imerode(morphedBW, strel('line', 5, 90));

    morphedBW = imclose(morphedBW, strel('rectangle', [5 12]));

    morphedBW = imerode(morphedBW, strel('line', 3, 0));
    morphedBW = imdilate(morphedBW, strel('line', 7, 90));
end

function morphedBW = morphBW(BW) % morphology operation for B W image
    morphedBW = imdilate(BW, strel('line', 9, 0));
    morphedBW = imerode(morphedBW, strel('line', 7, 90));

    morphedBW = imclose(morphedBW, strel('rectangle', [7 15]));

    morphedBW = imerode(morphedBW, strel('line', 5, 0));
    morphedBW = imdilate(morphedBW, strel('line', 5, 90));
end

% filter based on area and size
function idx = filterArea(props, imgSize)
    widthRatio = 0.85;
    heightRatio = 0.75;
    sizeRatio = 0.7;
    maxWidth = imgSize(2) * widthRatio;
    maxHeight = imgSize(1) * heightRatio;
    maxSize = (imgSize(1) * imgSize(2)) * sizeRatio;
    j=1;
    idx=0;
    for i=1:length(props)
        bb = props(i).BoundingBox;
        area=bb(3) * bb(4);
        if area> 500 && bb(3) > 30 && bb(4) > 10 && area < maxSize &&...
            bb(3) < maxWidth && bb(4) < maxHeight
            idx(j) = i;
            j=j+1;
        end
    end
end


function idx = filterAspectRatio(props) % check aspect ratio (w/h)
    j=1;
    idx=0;
    for i=1:length(props)
        bb = props(i).BoundingBox;
        aspectRatio=bb(3) / bb(4);
        if aspectRatio > 1.25 && aspectRatio < 10
            idx(j) = i;
            j=j+1;
        end
    end
end

% check centroid within selected region of image
function idx = filterPosition(props, imgSize) 
    imgHeight = imgSize(1);
    imgWidth = imgSize(2);

    minHeightRatio = 0.4;
    maxHeightRatio = 0.9;
    minWidthRatio = 0.25;
    maxWidthRatio = 0.75;
    xMin = imgWidth*(minWidthRatio);
    xMax = imgWidth*(maxWidthRatio);
    yMin = imgHeight*(minHeightRatio);
    yMax = imgHeight*(maxHeightRatio);
    xv = [xMin xMax xMax xMin xMin];
    yv = [yMin yMin yMax yMax yMin];

    hold;
    plot(xv,yv);

    j=1;
    idx=0;
    for i=1:length(props)
        centroid = props(i).Centroid;
        plot(centroid(1), centroid(2), 'r+');
        if inpolygon(centroid(1), centroid(2), xv, yv) > 0
            idx(j) = i;
            j=j+1;
        end
    end
end

function idx = filterCenter(props) % check centroid is within center area
    centerRatio = 0.35;
    minRatio = 0.5 - (centerRatio/2);
    maxRatio = 0.5 + (centerRatio/2);
    j=1;
    idx=0;
    for i=1:length(props)
        bb = props(i).BoundingBox;
        centroid = props(i).Centroid;
        xMin = bb(1) + bb(3)*(minRatio);
        xMax = bb(1) + bb(3)*(maxRatio);
        yMin = bb(2) + bb(4)*(minRatio);
        yMax = bb(2) + bb(4)*(maxRatio);
        centerxv = [xMin xMax xMax xMin xMin];
        centeryv = [yMin yMin yMax yMax yMin];
        if inpolygon(centroid(1), centroid(2), centerxv, centeryv) > 0
            idx(j) = i;
            j=j+1;
        end
    end
end

function filteredBW = findCCandFilter(BW) % find conn comp and filter
    cc = bwconncomp(BW);
    props = regionprops(cc, 'basic'); % area, centroid and bounding box

    if isempty(props)
        filteredBW = BW;
        return;
    end

    areaIdx = filterArea(props, size(BW));
    aspectRatioIdx = filterAspectRatio(props);
    centerIdx = filterCenter(props);
    positionIdx = filterPosition(props, size(BW));
    
    finalIdx = intersect(positionIdx, intersect(centerIdx, intersect(areaIdx, aspectRatioIdx)));
    filteredBW = ismember(labelmatrix(cc), finalIdx);
end

function finalcc = combineBRGY(BW1, BW2, BW3, BW4) % combine all cc
    cc1 = bwconncomp(BW1);
    cc2 = bwconncomp(BW2);
    cc3 = bwconncomp(BW3);
    cc4 = bwconncomp(BW4);

    finalcc = combineCC(cc1, cc2);
    finalcc = combineCC(finalcc, cc3);
    finalcc = combineCC(finalcc, cc4);
end

function cc = combineCC(cc1, cc2) % combine cc
    cc = cc1;
    cc.PixelIdxList = [cc.PixelIdxList [cc2.PixelIdxList]];
    cc.NumObjects = cc.NumObjects + cc2.NumObjects;
end

function averageOverlapRatio = checkOverlap(cc1, cc2) % calc overlapping ratio
    props1 = regionprops(cc1);
    props2 = regionprops(cc2);

    propsTable1 = struct2table(props1);
    propsBBoxArr1 = table2array(propsTable1(:, 3));

    propsTable2 = struct2table(props2);
    propsBBoxArr2 = table2array(propsTable2(:, 3));

    ratioUnion = bboxOverlapRatio(propsBBoxArr1, propsBBoxArr2);
    ratioMin = bboxOverlapRatio(propsBBoxArr1, propsBBoxArr2,"Min");

    averageOverlapRatio = (ratioUnion+ratioMin)/2;
end