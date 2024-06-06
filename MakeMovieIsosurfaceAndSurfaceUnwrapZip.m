function [] = MakeMovieIsosurfaceAndSurfaceUnwrapZip(data,CellID,tvec,savepath)
fs = 16;

clims = [10 100];

% Gaussian 3d filter
filtersize = [3 3 1];

% angle bins for cell surface
binsAngle = -180:5:180;
sigmaAngle = 2.5; % bin-center waiting sigma


zlayervec = 1:60; 
Nlayers = numel(zlayervec);


% structuring element used for dilations
SE1 = strel('disk',1);


pathImages = data.Source; % path to images
pathImages2 = data.Source_ArtifactSubtracted;
listM = data.ImageFileListMyo; % Image file list for Myosin


cd(pathImages)
cd ..
Tmatrix = load('gridAnalysis.mat').trackingMatrixZT;
CID1 = Tmatrix(:,7:8:end,5);
CID2 = Tmatrix(:,8:8:end,5);

zeromat = CID1 == 0;
CID1(zeromat) = NaN;
CID2(zeromat) = NaN;


CID1iscell = CID1 == CellID;
CID2iscell = CID2 == CellID;
cellrows = any(CID1iscell | CID2iscell, 2);
% get vertex positions to overlay on myosin surface heatmap
N1y = Tmatrix(cellrows,3:8:end,:);
N1x = Tmatrix(cellrows,4:8:end,:);
N2y = Tmatrix(cellrows,5:8:end,:)+N1y;
N2x = Tmatrix(cellrows,6:8:end,:)+N1x;


for t = tvec

    cd(pathImages2);
    VM = tiffreadVolume(listM{t});
    VM = imgaussfilt3(VM,filtersize);
    
    
    % get the Cell ID label-matrix to get cell surface positions
    cd(pathImages);
    cframefoldername = strcat('SegmentationData/',sprintf('frame%04d',t));   
    cd(cframefoldername);
    L = load('ImageBWlabel.mat').ImageBWlabel;
    
    % select z-layers specified by layervec.
    L = L(:,:,zlayervec);
    VM = VM(:,:,zlayervec);
    VM = double(VM);
    
    % if cell is not present (i.e. sum of logical cell vector ==0) skip
    % time point
    if sum(L==CellID,'all')==0
        continue;
    end


    SurfMatDataA = [];
    VertexMatDataA = [];
    for z=1:Nlayers
        
        IM = VM(:,:,z);
        
        cellmask = L(:,:,z)==CellID;
        % The cell mask is the region inside the watershed lines where
        % the membrane intensity is highest. To get the cell perimeter on
        % the watershed lines dilate the cell mask by one pixel and remove
        % the cellmask (or interior). Then we dilate the perimeter pixels
        % to get a thicker perimeter or better statistics for perimeter
        % intensity
        membranemask = imdilate(cellmask,SE1) & ~cellmask;
        membranemask = imdilate(membranemask,SE1);  
        
        if all(~cellmask(:)) % if the cell isn't present at all continue
            continue;
        end
        
        % for each layer find the centroid and the perimeter pixel locations
        S = regionprops(cellmask,'Centroid').Centroid;
        Cx = S(1);Cy = S(2);

        % find indices of membrane mask and find their angular coordinates
        ind = find(membranemask);
        N = numel(ind);
        [Ny,Nx] = ind2sub(size(IM),ind);
        Cxvec = Cx*ones(N,1);
        Cyvec = Cy*ones(N,1);
        Nx_c = Nx - Cxvec;
        Ny_c = Ny - Cyvec;
        AngleVec = atan2(Ny_c,Nx_c);
        AngleVecDeg_Surf = AngleVec*180/pi;
        Zvec_Surf = z*ones(N,1);
        IntVecM = IM(ind);
        SurfMatDataA = [SurfMatDataA;Zvec_Surf,AngleVecDeg_Surf, IntVecM];
        
        % Get Vertex data
        n1x = round(N1x(:,t,z));
        n1y = round(N1y(:,t,z));
        n2x = round(N2x(:,t,z));
        n2y = round(N2y(:,t,z));
        Nx = [n1x;n2x];
        Ny = [n1y;n2y];
        zerovec = Nx< 0.4;
        Nx(zerovec) = [];
        Ny(zerovec) = [];
        ind = sub2ind(size(IM),Ny,Nx);
        N = numel(Nx);
        cxvec = Cx*ones(N,1);
        cyvec = Cy*ones(N,1);
        Nx_c = Nx - cxvec;
        Ny_c = Ny - cyvec;
        AngleVec = atan2(Ny_c,Nx_c);
        AngleVecDeg_Vertex = AngleVec*180/pi;
        Zvec_Vertex = z*ones(N,1);
        IntVecM = IM(ind);
        VertexMatDataA = [VertexMatDataA;Zvec_Vertex,AngleVecDeg_Vertex, IntVecM];
        

    end
    

    
    close all

    subplot(1,2,1)
    IntensityIsosurface(L,VM,CellID);
    set(gca,'Position',[0.02 0.11 0.4 0.81]);
    caxis(clims)
    %title(strcat(num2str(round(t*spf/60)),' min'),'FontSize',fs)


    binsZ = zlayervec;
    sigmaZ = 0.5;

    subplot(1,2,2)
    [ result ] = interpolatePointVals2Mat_gauss( SurfMatDataA(:,[1 2 3]), binsZ, binsAngle, sigmaZ, sigmaAngle,0);
    imagesc(binsAngle,binsZ*0.2635,result);colorbar; %caxis(clims)
    ylabel('Depth (microns)','FontSize',fs)
    xlabel('Angle (degrees)','FontSize',fs)
    set(gca,'Xlim',[-180 180])
    xticks([-180 -90 0 90 180])
    %xticks([60 120])
    hold on
    scatter(VertexMatDataA(:,2),VertexMatDataA(:,1)*0.2635,9,'k','filled');hold off
    %set(gca,'Xlim',[45 135])
    colormap jet
    caxis(clims)
    set(gca,'Fontsize',20)
    %set(gca,'Position',[0.4 0.11 0.4 0.81]);
    pause
    
      
    %set(gcf,'color','w','Position',[167 380 842 417]);
    filename = strcat(savepath,'UnwrapCell',num2str(CellID),'Frame',num2str(t),'.tif'); 
    %print(gcf,filename,'-dtiff')
    %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end


    
end % function


function [] = IntensityIsosurface(L,ImMyosin,CellID)


L = flip(L,3);
ImMyosin = flip(ImMyosin,3);

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




% get the patch object using isosurface
p = isosurface(double(BWmembrane),0.5);
linearInd = sub2ind(size(BWmembrane), round(p.vertices(:,2)), round(p.vertices(:,1)),round(p.vertices(:,3)));

% select type of intensities to view: vertex or face intensities

VertexInts = ImMyosin(linearInd);
FaceInts = nanmean([VertexInts(p.faces(:,1)),VertexInts(p.faces(:,2)),VertexInts(p.faces(:,3))],2);



ColorVec = FaceInts;



% generate patch with properties
p2 = patch( ...
    'Vertices', p.vertices, ...
    'Faces', p.faces, ...
    'FaceVertexCData', ColorVec, ...
    'FaceColor', 'flat','EdgeColor','none');

% change things like aspect ratio, viewing angle, zoom
daspect([1 1 .10/.5])

az = -50;
el = 15; % elevation of 10 degrees above object 
view([az,el])
%zoom(1.05)
axis off
colormap jet
end

function [Row,Col] = FindCentroid(BW)

[Nrows,Ncols,Nlayers] = size(BW);
[ColGrid,RowGrid,~] = meshgrid(1:Ncols,1:Nrows,1:Nlayers);

Row = round(nanmean(RowGrid(BW)));
Col = round(nanmean(ColGrid(BW)));

end


function [CfitMat] = BleachingLimits(pathImages,listM,tvec,CIDS,zvec)
filtersize = [3 3 1];

tSampleVec = tvec(1):20:tvec(end);
N = numel(tSampleVec);
minVec = NaN(N,1);
maxVec = NaN(N,1);

ind = 1;
for t=tSampleVec
    cd(pathImages);
    VM = tiffreadVolume(listM{t});
    VM = double(VM(:,:,zvec));
    VM = imgaussfilt3(VM,filtersize);
    cframefoldername = strcat('SegmentationData/',sprintf('frame%04d',t));   
    cd(cframefoldername);
    L = load('ImageBWlabel.mat').ImageBWlabel(:,:,zvec);
    Mask = false(size(L));
    for ic=1:numel(CIDS)
        Mask = Mask | L==CIDS(ic) ;
    end
    Mask = bwperim(Mask);
    
    IntVec = VM(Mask);
    
    minVec(ind) = min(IntVec(:));
    maxVec(ind) = max(IntVec(:));
    ind = ind+1;
end   

fmin = fit(tSampleVec',minVec,'exp1');
fmax = fit(tSampleVec',maxVec,'exp1');

CminVec = feval(fmin,tvec);
CmaxVec = feval(fmax,tvec);


CfitMat = [CminVec,CmaxVec];


end % function
