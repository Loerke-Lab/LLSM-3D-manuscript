function [] = myosinCellImages_LS(data,cellNum,t,zlayer)

% this function is used to generate cell images of different depths corresponding to the 
% myosin cell isosurface and unwrapped map (should be boundary myosin only??!)


% move to directory with the data
cd(data.Source)
cd ..

% color for the interior of the cell
cellcolorRGB = [0 0 1];
%cellcolorRGB = [0.5 0 0.5];

% z-vector and z-step
zvec = 1:60;
zstep = 0.2635;
zvecMicrons = zvec*zstep;

for z = zlayer % this is the layer for the cells to be visualized at
    cd(data.Source)
    % imageGap = tiffreadVolume(data.ImageFileListGap43{t});
    % imageGap = im2double(imageGap(:,:,z));
    % E = strel('disk',15);
    % filtersize = 0.5;
    % IfiltGap = imgaussfilt3(imageGap,filtersize);
    %IfiltGap = imlocalbrighten(IfiltGap);
    %IfiltGap = imreducehaze(IfiltGap);

    imageGap = tiffreadVolume(data.ImageFileListGap43{t});
    imageGap = im2double(imageGap);
    filtersize = 0.5;
    IfiltGap = imgaussfilt3(imageGap,filtersize);
    IfiltGap = IfiltGap(:,:,z);


    % get myo images
    %cd(data.Source_ArtifactSubtracted)
    % imageMyo = tiffreadVolume(data.ImageFileListMyo{t});
    % imageMyo = im2double(imageMyo(:,:,z));
    % E = strel('disk',15);
    % filtersize = 0.5;
    % IfiltMyo = imgaussfilt3(imageMyo,filtersize);
    %IfiltMyo = imlocalbrighten(IfiltMyo);
    %IfiltMyo = imreducehaze(IfiltMyo);

    imageMyo = tiffreadVolume(data.ImageFileListMyo{t});
    imageMyo = im2double(imageMyo);
    filtersize = 0.5;
    IfiltMyo = imgaussfilt3(imageMyo,filtersize);
    IfiltMyo = IfiltMyo(:,:,z);


    [sx, sy] = size(IfiltGap);

    % get tracking label matrix
    cd(data.Source)
    cd('SegmentationData/')
    trackmatrixpath = strcat(sprintf('frame%04d',t),'/ImageBWlabel.mat');
    loadtrackmatrix = load(trackmatrixpath);
    Tmatrix = loadtrackmatrix.ImageBWlabel;
    TmatrixZ = Tmatrix(:,:,z);
    cellmask = TmatrixZ == cellNum;

    RGBcells1 = zeros(sx,sy);
    RGBcells2 = zeros(sx,sy);
    RGBcells3 = zeros(sx,sy);
    % create color image of cells
    BWtrackcells = zeros(sx,sy);

    RGBcells1(cellmask) = cellcolorRGB(1);
    RGBcells2(cellmask) = cellcolorRGB(2);
    RGBcells3(cellmask) = cellcolorRGB(3);

    BWtrackcells(cellmask) = 1;

    RGBcells = cat(3,RGBcells1,RGBcells2,RGBcells3);

    %image1Gap = imageGap;
     image1Gap = IfiltGap;
%     image1Gap = imtophat(image1Gap,E);
    %image1Gap = imlocalbrighten(image1Gap);
    % minGap = min(image1Gap,[],'all');
    % maxGap = max(image1Gap,[],'all');
    %image1Gap = imadjust(image1Gap);

    image3Gap = repmat(image1Gap,[1,1,3]);

    %image1Myo = imageMyo;
     image1Myo = IfiltMyo;
%     image1Myo = imtophat(image1Myo,E);

    %image1Myo = imlocalbrighten(image1Myo,0.1);
    % minMyo = min(image1Myo,[],'all');
    % maxMyo = max(image1Myo,[],'all');
    %image1Myo = imadjust(image1Myo);
    

    image3Myo = repmat(image1Myo,[1,1,3]);

%     rows = size(image1Myo,1);
%     cols = size(image1Myo,2);
%     image3Myo(:,:,1) = image1Myo;
%     image3Myo(:,:,2) = zeros(rows,cols);
%     image3Myo(:,:,3) = zeros(rows,cols);;

    % transpose the color segmented image on the grayscale image
    imageRGBGap = 0.75*image1Gap + 0.25*RGBcells;

    % find centroid of the collection of the cells associated with the
    % interfaces. First, get BW of cells.
    stats = regionprops(BWtrackcells,'Centroid','BoundingBox');
    Centroid = stats.Centroid;

    % bounding box gives the smallest rectangle that holds all the cells. Since
    % we would like a square window we use the larger of the two dimensions and
    % enlarge it by 50%
    BoundingBox = stats.BoundingBox;

    % update cropsize  
    %if t == tvec(1)
        cropsize = 200;%round(2*max(BoundingBox(3:4)));
    %end

    % crop image
    rect = [round(Centroid(1) - cropsize/2) round(Centroid(2) - cropsize/2) cropsize cropsize];
    IcropGap = imcrop(imageRGBGap,rect);
    IcropMyo = imcrop(image1Myo,rect);
    %Icrop = imrotate(Icrop,90);

    %IcropGap = imadjust(IcropGap);
    IcropMyo = imadjust(IcropMyo);

    RGBCellsCrop = imcrop(imageRGBGap,rect);
    image1GapCrop = imcrop(image1Gap,rect);
    image1GapCrop = imadjust(image1GapCrop);

    IcropGap = 0.75*image1GapCrop + 0.25*RGBCellsCrop;

    %imshow(IcropGap);
    figure
    subplot(1,2,1)
    imshow(IcropGap)
    subplot(1,2,2)
    imshow(IcropMyo)


end



