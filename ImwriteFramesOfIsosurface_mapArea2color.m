function [] = ImwriteFramesOfIsosurface_mapArea2color(data,cellrownum,frameSelection,savepath)
% 

% if just wanting to view images but not create a movie, set to false.
pathImages = data.Source;
zvec = 1:60;
spf = data.SecPerFrame; 
XYsf = 0.1;
Zsf = 0.5;


% Load Geometry data
cd(pathImages)
cd ..
loadGdata = load('GeometryData.mat');
CellNumbers = loadGdata.CellNumbers;
CellID = CellNumbers(cellrownum);
CellAreaZvT = squeeze(loadGdata.Area(cellrownum,:,zvec))'*(XYsf)^2;
frameStart = find(all(isfinite(CellAreaZvT),1),1,'first');
frameEnd = find(all(isfinite(CellAreaZvT),1),1,'last');
frameVec = frameStart:frameEnd;
tvecMinArea = frameVec*spf/60;
CellAreaZvT = CellAreaZvT(:,frameVec);




[AreaFilt] = filterImage3DpaddedEdges(CellAreaZvT, 'Gauss', [2,1,0]);


jetmap = jet(256);
imagesc(AreaFilt);colorbar;colormap(jetmap)
pause;

% maxR = max(AreaFilt(:));
% minR = min(AreaFilt(:));
% clims = [minR maxR];


AreaFilt = AreaFilt(:,frameSelection);
maxR = max(AreaFilt(:));
minR = min(AreaFilt(:));
clims = [minR maxR];

frameVec= frameVec(frameSelection);
clf;

% If you want to change the view angle use linspace, otherwise set to
% constant
FrameNum = numel(frameVec);
frameIdx = 1;
for it = 1:FrameNum
    
    t = frameVec(it);
    RateAtTimePointVec = AreaFilt(:,it);
    RateAtTimePointVec = reshape(RateAtTimePointVec,[1 1 numel(RateAtTimePointVec)]);
    
    
    
    % load Current image stack
    cd(pathImages);
    cframefoldername = strcat('SegmentationData/',sprintf('frame%04d',t));   
    cd(cframefoldername);
    L = load('ImageBWlabel.mat').ImageBWlabel;
    L = L(:,:,zvec);
    [sx,sy,~] = size(L);
    
    IntensityMat = repmat(RateAtTimePointVec,[sx,sy,1]);
    

    L = flip(L,3);
    IntensityMat = flip(IntensityMat,3);
%     
%     subplot(1,2,1)
    IsosurfaceWithIntensity(L,IntensityMat,CellID,jetmap,clims);
   % set(gca,'Position',[0.05    0.1100    0.3347    0.8150])
%     
%     subplot(1,2,2);
%     imagesc(tvecMinRate,zvecMicrons,AreaRateFilt)
%     hold on
%     plot([t*spf/60 t*spf/60],[1 30],'k--')
%     xlabel('Time(min)','FontSize',fs,'Color','k')
%     ylabel('Depth (microns)','FontSize',fs,'Color','k')
%     title('Cell Area Rate (microns^2/min)','FontSize',fs,'Color','k')
%     colormap(b2rmap)
%     caxis(clims);
%     c= colorbar;
%     %c.Color = 'w';
%     set(gca,'Position',[0.403    0.1134    0.5    0.7997])
%     set(gca,'Color','w')
%     set(gcf,'Position',[260   478   802   319])
     view(0,10)
     if nargin > 3
        frame = getframe(gcf);   
        filename = strcat(savepath,'Frame',num2str(frameIdx),'.tif');
        RGB = frame2im(frame);
        imwrite(RGB,filename);
        frameIdx = frameIdx + 1;
     else
         pause;
    end
    clf


end


end


function [] = IsosurfaceWithIntensity(L,ImMyosin,CellID,colorlist,clims)
% select the cell of interest
BWcell = L == CellID;


% since the segmentation excludes the boundary lines dilate by one pixel to get membrane 
se1 = strel('disk',1);
BWmembrane = imdilate(BWcell,se1);


% select out the centroid of the cell, then crop a region around the cell
[Row,Col] = FindCentroid(BWcell);
BoxHalfSize = 100;
rowvec = (Row-BoxHalfSize):(Row+BoxHalfSize);
colvec = (Col-BoxHalfSize):(Col+BoxHalfSize);
BWmembrane = BWmembrane(rowvec,colvec,:);
ImMyosin= ImMyosin(rowvec,colvec,:);


ImMyosin= im2double(ImMyosin);
filtersize = [5 5 1];
ImMyosin = imgaussfilt3(ImMyosin,filtersize);

% get the patch object using isosurface
p = isosurface(double(BWmembrane),0.5);
linearInd = sub2ind(size(BWmembrane), round(p.vertices(:,2)), round(p.vertices(:,1)),round(p.vertices(:,3)));

% select type of intensities to view: vertex or face intensities

VertexInts2 = ImMyosin(linearInd);
FaceInts2 = nanmean([VertexInts2(p.faces(:,1)),VertexInts2(p.faces(:,2)),VertexInts2(p.faces(:,3))],2);



ColorVec2 = FaceInts2;


figure(1)
% generate patch with properties
p2 = patch( ...
    'Vertices', p.vertices, ...
    'Faces', p.faces, ...
    'FaceVertexCData', ColorVec2, ...
    'FaceColor', 'flat','EdgeColor','none','FaceAlpha',0.6);

% change things like aspect ratio, viewing angle, zoom
daspect([1 1 .10/.5])

az = -50;
el = 15; % elevation of 10 degrees above object 
view([az,el])
axis off
colormap(colorlist)
caxis(clims)
colorbar 
set(gcf,'Color','w')


end

function [Row,Col] = FindCentroid(BW)

[Nrows,Ncols,Nlayers] = size(BW);
[ColGrid,RowGrid,~] = meshgrid(1:Ncols,1:Nrows,1:Nlayers);

Row = round(nanmean(RowGrid(BW)));
Col = round(nanmean(ColGrid(BW)));

end

function [Rate,tvec_centered] = getRate(Mat,DeltaT)
Nframes = size(Mat,2);
tvec2 = (DeltaT+1):Nframes;
tvec1 = 1:(Nframes-DeltaT);
Mat2 = Mat(:,tvec2);
Mat1 = Mat(:,tvec1);
Rate = (Mat2 - Mat1)/DeltaT;
tvec_centered = (tvec1+tvec2)/2;
end

