function [] = SingleCell3DIsosurfaceMovie_withCrossSections(data,MovieNum,CellNumber,tvec)
% This function generates a 3D isosurface for a single cell.
% This function was used to produce the single-cell isosurface and
% cross-sections for the Figure 1. While a time vector can be an input,
% each time point overwrites the previous time point. I used it for a
% single time point as in:
%  SingleCell3DIsosurfaceMovie_withCrossSections(data,1,3,1)
%
% 4/29/22


close all

% structuring element for erosion of cell mask, for the purpose of coloring
% the cell inside and not along the plasma membrane
SE = strel('disk',4);

% vector of layers to print cross-section images of
layervec = linspace(1,85,7);


% list of image file paths
list = data(MovieNum).ImageFileListGap43;

pathstr = data(MovieNum).Source;

% Move to directory
cd(data(MovieNum).Source)
cd ..


% color to color the cell
ColorRGB = [0 1 1]; % greenish-blue


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
    
    % generate a 3D isosurface and print eps to desktop
    epithelialImage3D(LabelImage,CellNumber,ColorRGB,1);
    print(gcf,'~/Desktop/CellIsosurface.eps','-depsc')
     
    
    % print out cross-sectional images for layers specified by layervec
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


        cellmask = Tmatrix(:,:,layer) == CellNumber;
        if all(~cellmask,'all')
            continue;
        end
        cellmask = imerode(cellmask,SE);

        RGBcells1(cellmask) = ColorRGB(1);
        RGBcells2(cellmask) = ColorRGB(2);
        RGBcells3(cellmask) = ColorRGB(3);

        
        % find centroid of the collection of the cells associated with the
        % interfaces. First, get BW of cells.
        statsA = regionprops(cellmask,'Centroid','BoundingBox');
        CentroidA = statsA.Centroid;
        cropsize = 140;
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
        filename = strcat('~/Desktop/SampleVideo/Layer',num2str(layer),'.tif');
        imwrite(imageRGB,filename);
    end

    
 
    
end % time loop   
end % function


function [] = epithelialImage3D(LabelImage,CellNumber,colorsArray,smoothflag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all


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



