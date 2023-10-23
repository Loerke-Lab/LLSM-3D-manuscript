function [] = ImwriteFramesOfIsosurface_mapAreaRate2color_Rho_Wire(data,cellrownum,frameSelection,savepath)

% if just wanting to view images but not create a movie, set to false.
pathImages = data.Source;
zvec = 4:19;
spf = data.SecPerFrame; 
DeltaT = round(30/spf);
XYsf = 0.163;



% Load Geometry data
cd(pathImages)
%cd ..
loadGdata = load('GeometryData.mat');
CellNumbers = loadGdata.CellNumbers;
CellID = CellNumbers(cellrownum);
CellAreaZvT = squeeze(loadGdata.Area(cellrownum,:,zvec))'*(XYsf)^2;
frameStart = find(all(isfinite(CellAreaZvT),1),1,'first');
frameEnd = find(all(isfinite(CellAreaZvT),1),1,'last');
frameVec = frameStart:frameEnd;
CellAreaZvT = CellAreaZvT(:,frameVec);


CellAreaZvT = filterImage3DpaddedEdges(CellAreaZvT, 'Gauss', 2);


% Get the rate of change
[dAdtmat,framevec_centered] = getRate(CellAreaZvT,frameVec,DeltaT);

dAdtmat = dAdtmat/spf*60; % convert rate to square microns per minute.
framevec_centered = round(framevec_centered); % round because the centered vec can be non-integer.



maxR = max(dAdtmat(:));
minR = abs(min(dAdtmat(:)));
maxmaxR = max(maxR,minR);
clims = [-maxmaxR maxmaxR];
%clims = [-2 2];
b2r_cmap = b2r(-maxmaxR,maxmaxR);
b2r_cmap = flip(b2r_cmap,1);
imagesc(dAdtmat);colorbar;colormap(b2r_cmap)
pause;
%findices =

dAdtmat = dAdtmat(:,frameSelection);
framevec_centered = framevec_centered(frameSelection);


maxR = max(dAdtmat(:));
minR = abs(min(dAdtmat(:)));
maxmaxR = max(maxR,minR);
clims = [-maxmaxR maxmaxR];
%clims = [-2 2];
b2r_cmap = b2r(-maxmaxR,maxmaxR);
b2r_cmap = flip(b2r_cmap,1);
clf;


% If you want to change the view angle use linspace, otherwise set to
% constant
FrameNum = numel(frameSelection);
frameIdx = 1;
for it = 1:FrameNum
    
    t = framevec_centered(it);
    %% 
    RateAtTimePointVec = dAdtmat(:,it);
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
    IsosurfaceWithIntensity(L,IntensityMat,CellID,b2r_cmap,clims);
    %view(-83,10)
    view(160,10)
    
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
    
     if nargin > 3
        frame = getframe(gcf);   
        filename = strcat(savepath,'Frame',num2str(frameIdx),'.tif');
        %filename = strcat(savepath,'Frame',num2str(frameIdx),'.eps');
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
%BoxHalfSize = 100;
BoxHalfSize = 75; % added KB
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
%caxis(clims)
caxis([-6 6])
colorbar 
set(gcf,'Color','w')


end

function [Row,Col] = FindCentroid(BW)

[Nrows,Ncols,Nlayers] = size(BW);
[ColGrid,RowGrid,~] = meshgrid(1:Ncols,1:Nrows,1:Nlayers);

Row = round(nanmean(RowGrid(BW)));
Col = round(nanmean(ColGrid(BW)));

end

function [Rate,fvec_centered] = getRate(Mat,frameVec,DeltaT)
Nframes = size(Mat,2);
tvec2 = (DeltaT+1):Nframes;
tvec1 = 1:(Nframes-DeltaT);
Mat2 = Mat(:,tvec2);
Mat1 = Mat(:,tvec1);
Rate = (Mat2 - Mat1)/DeltaT;

fvec2 = frameVec((DeltaT+1):Nframes);
fvec1 = frameVec(1:(Nframes-DeltaT));
fvec_centered = (fvec1+fvec2)/2;
end

