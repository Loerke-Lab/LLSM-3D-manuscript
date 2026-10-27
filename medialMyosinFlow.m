function [] = medialMyosinFlow(data,MovieNum,CellNumber,tVec)

% this function produces example cell images of medial myosin flows toward vertices


% structuring element for erosion of cell mask, for the purpose of coloring
% the cell inside and not along the plasma membrane
SE = strel('disk',4);

zsf = 0.2635;

layer = round(5/zsf);
%layer = 1;

% list of image file paths
listGap43 = data(MovieNum).ImageFileListGap43;
listMyo = data(MovieNum).ImageFileListMyo;

pathstr = data(MovieNum).Source;

% Move to directory
cd(data(MovieNum).Source)
cd ..


% color to color the cell
ColorRGB = [0 1 1]; % greenish-blue

% file name and path for saving the movie. Note: Here I save it to my
% desktop.
filename = strcat('~/Desktop/Movie_',num2str(MovieNum),'_Cell_',num2str(CellNumber),'_5um.mp4');
writerObj = VideoWriter(filename,'MPEG-4');
%writerObj.FrameRate = 15; % This can be changed to alter the playback speed.
writerObj.FrameRate = 7.5; % added 4/26/22 KB
writerObj.Quality = 100;
open(writerObj);

for t=tVec


    % get tracking label matrix
    trackmatrixpath = strcat(pathstr,'/SegmentationData/',sprintf('frame%04d',t),'/ImageBWlabel.mat');
    loadtrackmatrix = load(trackmatrixpath);
    Tmatrix = loadtrackmatrix.ImageBWlabel;
    LabelImage = Tmatrix;


    [sx,sy,~] = size(LabelImage);
    zeromat = zeros(sx,sy);
    LabelImage = cat(3,zeromat,LabelImage,zeromat);
    LabelImage = flip(LabelImage,3);


    cd(pathstr)
    E = strel('disk',20);


    % get cell membrane image
    V = tiffreadVolume(listGap43{t});
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
    %rect = [round(CentroidA(1) - cropsize/2) round(CentroidA(2) - cropsize/2) cropsize cropsize];
    rect = [floor(CentroidA(1) - cropsize/2) floor(CentroidA(2) - cropsize/2) cropsize cropsize];
    RGBcells1 = imcrop(RGBcells1,rect);
    RGBcells2 = imcrop(RGBcells2,rect);
    RGBcells3 = imcrop(RGBcells3,rect);
    image1 = imcrop(image1,rect);

    RGBcells = cat(3,RGBcells1',RGBcells2',RGBcells3');

    
    image1 = imtophat(image1,E)';
    image1 = imadjust(image1);

    image3 = repmat(image1,[1,1,3]);
    imageRGB = 0.60*image3 + 0.40*RGBcells;


    % get myosin image
    V = tiffreadVolume(listMyo{t});
    V = im2double(V);
    filtersize = 0.5;
    Ifilt = imgaussfilt3(V,filtersize);
    image2 = Ifilt(:,:,layer);

    image2 = imcrop(image2,rect);

    image2 = imtophat(image2,E)';
    image2 = imadjust(image2);


    figure(1)
    % subplot(1,2,1)
    % imshow(image1)
    % subplot(1,2,2)
    % imshow(image2)

    %Image3D2(:,:,1) = image2;
    %Image3D2(:,:,2) = zeros(size(image2,1),size(image2,2));
    %Image3D2(:,:,3) = zeros(size(image2,1),size(image2,2));

    image12 = imfuse(image1,image2);%,'Scaling','none');%,'ColorChannels',[2 1 0]);
    imshow(image12)
    %imshow(Image3D2)
    %axis([round(CentroidA(1) - cropsize/2) round(CentroidA(1) + cropsize/2) round(CentroidA(2) - cropsize/2) round(CentroidA(2) + cropsize/2)])

    %pause(0.01)

    frame = getframe(gcf);
    writeVideo(writerObj,frame);



end
