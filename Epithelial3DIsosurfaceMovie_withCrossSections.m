function [] = Epithelial3DIsosurfaceMovie_withCrossSections(data,MovieNum)
% This function generates a 3D isosurface of the entire tissue with random
% coloring of the cells and also 3 cross-sectional images with coloring of
% the cell interiors
%
% Example:
% Epithelial3DIsosurfaceMovie_withCrossSections(dataLSTC,3).


% 4/26/22 Tim Vanderleest

close all
% structuring element to erode cell masks to color the cell interiors
SE = strel('disk',4);


% layers to show produce cross-sections of
ApicalLayer = 1;
LateralLayer = 30;
BasalLayer = 60;


% list of image file paths
list = data(MovieNum).ImageFileListGap43;

pathstr = data(MovieNum).Source;

% Move to directory
cd(data(MovieNum).Source)
cd ..

% get full cell label Vector
loadGeometry = load('GeometryData.mat');
CellNumbers = loadGeometry.CellNumbers;
NumCells = length(CellNumbers);

% Generate colors for the cells (Nx3 array for N cells).
% Randomly permute the rows in the colors array so that cells have random
% colors instead of the tissue having a color gradient.
ColorsArray = jet(NumCells);
p = randperm(NumCells);
ColorsArray(p,:) = ColorsArray;


%time vector, have it do the entire movie or the size of the list
tvec = 1:size(list,1);
frameIdx = 1;
for t = tvec

    
    % get tracking label matrix
    trackmatrixpath = strcat(pathstr,'/SegmentationData/',sprintf('frame%04d',t),'/ImageBWlabel.mat');
    loadtrackmatrix = load(trackmatrixpath);
    Tmatrix = loadtrackmatrix.ImageBWlabel;
    LabelImage = Tmatrix;
    
    %[LabelDomes] = createCellDomes(LabelImage(:,:,1));
    
    [sx,sy,~] = size(LabelImage);
    zeromat = zeros(sx,sy);
    
    % pad the top and bottom of the matrix with zeros, otherwise we get
    % holes in the top/bottom.
    LabelImage = cat(3,zeromat,LabelImage,zeromat);
    LabelImage = flip(LabelImage,3);
    

    epithelialImage3D(LabelImage,CellNumbers,ColorsArray,1);
    frame = getframe(gcf);   
    filename = strcat('~/Desktop/Isosurface',num2str(frameIdx),'.tif');
    RGB = frame2im(frame);
    imwrite(RGB,filename);
    
    

    cd(pathstr)
    V = tiffreadVolume(list{t});
    V = im2double(V);
    E = strel('disk',15);
    filtersize = 0.5;
    Ifilt = imgaussfilt3(V,filtersize);
    
    %% Apical Layer cross-section
    % create color image of cells
    RGBcells1 = zeros(sx,sy);
    RGBcells2 = zeros(sx,sy);
    RGBcells3 = zeros(sx,sy);
    for k=1:NumCells
        cellmask = Tmatrix(:,:,ApicalLayer) == CellNumbers(k);
        if all(~cellmask,'all')
            continue;
        end
        cellmask = imerode(cellmask,SE);
        
        RGBcells1(cellmask) = ColorsArray(k,1);
        RGBcells2(cellmask) = ColorsArray(k,2);
        RGBcells3(cellmask) = ColorsArray(k,3);
    end
    RGBcells = cat(3,RGBcells1',RGBcells2',RGBcells3');
    
    image1 = Ifilt(:,:,ApicalLayer)';
    image1 = imtophat(image1,E);
    image1 = imadjust(image1);
    
    image3 = repmat(image1,[1,1,3]);
    imageRGB = 0.60*image3 + 0.40*RGBcells;    
    imwrite(imageRGB,'~/Desktop/ApicalSlice.tif')
    
    
    %% Lateral Layer cross-section
    % create color image of cells
    RGBcells1 = zeros(sx,sy);
    RGBcells2 = zeros(sx,sy);
    RGBcells3 = zeros(sx,sy);
    for k=1:NumCells
        cellmask = Tmatrix(:,:,LateralLayer) == CellNumbers(k);
        if all(~cellmask,'all')
            continue;
        end
        cellmask = imerode(cellmask,SE);
        
        RGBcells1(cellmask) = ColorsArray(k,1);
        RGBcells2(cellmask) = ColorsArray(k,2);
        RGBcells3(cellmask) = ColorsArray(k,3);
    end
    RGBcells = cat(3,RGBcells1',RGBcells2',RGBcells3');
    
    image1 = Ifilt(:,:,LateralLayer)';
    image1 = imtophat(image1,E);
    image1 = imadjust(image1);
    
    image3 = repmat(image1,[1,1,3]);
    imageRGB = 0.60*image3 + 0.40*RGBcells;
    imwrite(imageRGB,'~/Desktop/LateralSlice.tif')
    
    
    
    %% Basal Layer cross-section
    % create color image of cells
    RGBcells1 = zeros(sx,sy);
    RGBcells2 = zeros(sx,sy);
    RGBcells3 = zeros(sx,sy);
    for k=1:NumCells
        cellmask = Tmatrix(:,:,BasalLayer) == CellNumbers(k);
        if all(~cellmask,'all')
            continue;
        end
        cellmask = imerode(cellmask,SE);
        
        RGBcells1(cellmask) = ColorsArray(k,1);
        RGBcells2(cellmask) = ColorsArray(k,2);
        RGBcells3(cellmask) = ColorsArray(k,3);
    end
    RGBcells = cat(3,RGBcells1',RGBcells2',RGBcells3');
    
    image1 = Ifilt(:,:,BasalLayer)';
    image1 = imtophat(image1,E);
    image1 = imadjust(image1);
    
    image3 = repmat(image1,[1,1,3]);
    imageRGB = 0.60*image3 + 0.40*RGBcells;
    imwrite(imageRGB,'~/Desktop/BasalSlice.tif')
    
    pause;

    
end % time loop   
end % function


function [] = epithelialImage3D(LabelImage,CellNumbersVec,colorsArray,smoothflag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all


% set all negative indices (if they exist) to zeros
negativeIndices = LabelImage == -1;
LabelImage(negativeIndices) = 0;

% get cell labels and number of cells
PositiveLabels = LabelImage > 1;
allLabels = unique(LabelImage(PositiveLabels));


NumCells = length(allLabels);

% flip image upside-down so that apical side is up.
%LabelImage = flipdim(LabelImage,3);

for jj=1:NumCells
    label = allLabels(jj);
    
    if sum(LabelImage(:)== label) == 0
        continue
    end
    % 3D data smoothing
    if smoothflag
        p = patch(isosurface(smooth3(LabelImage==label),0.5));
        %p = patch(isosurface(smooth3(LabelImage==label,'gaussian',[9 9 5],3),0.5));
    else
        p = patch(isosurface(LabelImage==label));
    end
    daspect([0.263,0.263,0.10])     

    set(p,'FaceColor',colorsArray(ismember(CellNumbersVec,label),:),'EdgeColor','none','FaceAlpha',0.6);
    

    
    hold on
end
view(90,20); %axis([1 158 1 397 1 52])
%light
camlight 
%lighting gouraud

axis off
set(gcf,'color','w');
text


% % find center of mass of all cells and center axis on that point
% BWallcells = zeros(size(LabelImage));
% for i=1:length(allLabels)
%     BWallcells(LabelImage==allLabels(i)) = 1;
% end
% STATS = regionprops(BWallcells,'Centroid');
% Cent = STATS.Centroid;
% axis([Cent(1)-200 Cent(1)+200 Cent(2)-200 Cent(2)+200 4 10])
% hold off
% axis off


end


function [BWdome] = createCellDomes(LabelImage2D)
maxLayers = 5;
BWdome = zeros([size(LabelImage2D),maxLayers]);


BW = LabelImage2D == 0;

Dist = bwdist(BW);

DistSQRT = sqrt(Dist);

for ii=1:maxLayers
    maskIsoDist = DistSQRT >= ii;
    labelIsoDist = LabelImage2D;
    labelIsoDist(~maskIsoDist) = 0;
    
    BWdome(:,:,maxLayers-ii+1) = labelIsoDist;
end    
BWdome(:,:,1) = [];
BWdome(:,:,end) = [];
end



