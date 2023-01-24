function [] = FourCell3DIsosurfaceMovie_withCrossSections(data,MovieNum,CellNumbers,tvec)
% This function generates a 3D isosurface for a 4-cell group (i.e. T1 group).
% This function was used to produce a 3D isosurface and cross-sections. 
% While a time vector can be an input, each time point overwrites the 
% previous time point. I used it for a single time point as in:
%  FourCell3DIsosurfaceMovie_withCrossSections(dataSample,1,[3 4 8 9],1)
%
% 6/20/22

close all
SE = strel('disk',4);
layervec = linspace(1,85,7);


% list of image file paths
list = data(MovieNum).ImageFileListGap43;

pathstr = data(MovieNum).Source;

% Move to directory
cd(data(MovieNum).Source)
cd ..

% get full cell label Vector
%loadGeometry = load('GeometryData.mat');
%CellNumbers = loadGeometry.CellNumbers;


NumCells = 1;
%ColorRGB = [0 1 1]; % greenish-blue
ColorRGBs = [0 1 1;1 0 0;0 0 1;0 1 0];

% Generate colors for the cells (Nx3 array for N cells).
% Randomly permute the rows in the colors array so that cells have random
% colors instead of the tissue having a color gradient.
% ColorsArray = jet(NumCells);
% p = randperm(NumCells);
% ColorsArray(p,:) = ColorsArray;
%ColorsArray = [0 0 1];


%time vector
for t = tvec

    
    % get tracking label matrix
    trackmatrixpath = strcat(pathstr,'/SegmentationData/',sprintf('frame%04d',t),'/ImageBWlabel.mat');
    loadtrackmatrix = load(trackmatrixpath);
    Tmatrix = loadtrackmatrix.ImageBWlabel;
    LabelImage = Tmatrix;
    
    
    [sx,sy,~] = size(LabelImage);
    zeromat = zeros(sx,sy);
    LabelImage = cat(3,zeromat,LabelImage,zeromat);
    LabelImage = flip(LabelImage,3);
    
    for ii=1:4
        epithelialImage3D(LabelImage,CellNumbers(ii),ColorRGBs(ii,:),1);hold on
    end
    camlight
    print(gcf,'~/Desktop/CellIsosurface.eps','-depsc')
%     frame = getframe(gcf);   
%     filename = strcat('~/Desktop/3D Figures/Isosurface',num2str(frameIdx),'.tif');
%     RGB = frame2im(frame);
%     imwrite(RGB,filename);
%     
    
    
    cd(pathstr)
    E = strel('disk',20);
    for layer = layervec

        
        V = tiffreadVolume(list{t});
        V = im2double(V);
        filtersize = 0.5;
        Ifilt = imgaussfilt3(V,filtersize);
        image1 = Ifilt(:,:,layer);

        
        
        % create color image of cells
        RGBcells1 = zeros(sx,sy);
        RGBcells2 = zeros(sx,sy);
        RGBcells3 = zeros(sx,sy);
        
        cellmaskAll = false(sx,sy);
        for ii=1:4
            cellmask = Tmatrix(:,:,layer) == CellNumbers(ii);
            cellmaskAll = cellmaskAll | cellmask;
            if all(~cellmask,'all')
                continue;
            end
            cellmask = imerode(cellmask,SE);

            RGBcells1(cellmask) = ColorRGBs(ii,1);
            RGBcells2(cellmask) = ColorRGBs(ii,2);
            RGBcells3(cellmask) = ColorRGBs(ii,3);
        end

        
        % find centroid of the collection of the cells associated with the
        % interfaces. First, get BW of cells.
        statsA = regionprops(cellmaskAll,'Centroid','BoundingBox');
        CentroidA = NaN(4,2);
        for ii=1:4
            CentroidA(ii,:) = statsA(ii).Centroid;
        end
        CentroidA = mean(CentroidA);
        cropsize = 200;
        % crop image
        rect = [round(CentroidA(1) - cropsize/2) round(CentroidA(2) - cropsize/2) cropsize cropsize];
        RGBcells1 = imcrop(RGBcells1,rect);
        RGBcells2 = imcrop(RGBcells2,rect);
        RGBcells3 = imcrop(RGBcells3,rect);
        image1 = imcrop(image1,rect);
        
        
        
        RGBcells = cat(3,RGBcells1',RGBcells2',RGBcells3');

        
        image1 = imtophat(image1,E)';
        image1 = imadjust(image1);

        image3 = repmat(image1,[1,1,3]);
        imageRGB = 0.60*image3 + 0.40*RGBcells;
        filename = strcat('~/Desktop/Layer',num2str(layer),'.tif');
        imwrite(imageRGB,filename);
    end

    
 
    
end % time loop   
end % function


function [] = epithelialImage3D(LabelImage,CellNumber,colorsArray,smoothflag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% set all negative indices (if they exist) to zeros
negativeIndices = LabelImage == -1;
LabelImage(negativeIndices) = 0;

% get cell labels and number of cells
PositiveLabels = LabelImage > 1;
allLabels = unique(LabelImage(PositiveLabels));


NumCells = length(CellNumber);

% flip image upside-down so that apical side is up.
%LabelImage = flipdim(LabelImage,3);

for jj=1:NumCells
    label = CellNumber(jj);
    
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

    set(p,'FaceColor',colorsArray,'EdgeColor','none','FaceAlpha',0.6);
    

    
    hold on
end
view(90,15); %axis([1 158 1 397 1 52])
%light
%camlight 
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



