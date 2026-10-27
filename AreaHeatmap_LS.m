function [] = AreaHeatmap_LS(data,MovieNums,startCell)
% AreaHeatmap_LS  Produces a heatmap of the cell area
%
% Input:       data: structure that contains the following fields:                  
%                   -'ImageFileList': a list of file names for the
%                     cropped images.
%                   -'Source': the path to the cropped images.
%                   -'SecPerFrame': number of seconds per time frame
%
%               MovieNums: integer vector of which movies to process.
% 
%               startCell: integer scalar of cell index to start the loop iteration with 
%
%               
%
% Output:       No outputs, just displays a figure and pauses
%

od = cd; % get original directory

fs = 16; % font size for figures

% length scale factor for X-Y (microns per pixel)
lsf = 0.104;

% length scale factor for Z (mircons per pixel)
zsf = 0.2635;
zvec = 1:70; % z-layers to select from z-stack
zvecMic = zvec*zsf;

Nmovies = numel(MovieNums);
for mn = 1:Nmovies
    
    % Move to directory and load required data
    cd(data(MovieNums(mn)).Source)
    cd ..

    spf = data.SecPerFrame; % seconds per frame conversion factor
    
    % Load cell Area matrix and convert to microns^2
    loadGdata = load('GeometryData.mat');
    AreaDataMat = loadGdata.Area(:,:,zvec)*lsf^2;

    NumCells = size(AreaDataMat,1);

    for int=startCell:NumCells

        CellArea = squeeze(AreaDataMat(int,:,:))';
        finite = isfinite(CellArea);
        startFrame = find(all(finite,1),1,'first');
        stopFrame   = find(all(finite,1),1,'last');

        %tvec = startFrame:stopFrame;

        [~, tvec] = FiniteSequenceFromNans(mean(CellArea,1));%,1);
        
        if isempty(tvec)
            continue;
        end


        AreaA = CellArea(:,tvec);
      
        % smooth out the area matrix
        AreaAF = filterImage3DpaddedEdges(AreaA, 'Gauss', 2);
        %AreaAF = filterImage3DpaddedEdges(AreaA, 'Gauss', [2,1,0]);

        tvecMin = (1:size(AreaAF,2))*spf/60;
        
        if isempty(tvecMin) || tvecMin(end) < 5
            continue;
        end
        imagesc(tvecMin,zvecMic,AreaAF);
        colormap(parula);
        colorbar
        ylabel('Z-depth (microns)','FontSize',fs)
        xlabel('Time (min)','FontSize',fs)
        set(gca,'FontSize',fs)
        title(num2str(int),'FontSize',fs)
        maxL = max(abs(AreaAF),[],'all');
        minL = min(abs(AreaAF),[],'all');
        clim(gca,[minL, maxL]);
        
        
        choice = menu('Choose one: ','Select time window','Skip');
        if choice ==1
            StartStop = input('Input start and stop vec in min: ');
            
            
            idx = tvecMin>=StartStop(1) & tvecMin <= StartStop(2);
            
            
            AreaAF = AreaAF(:,idx);

            
            tvecMin = (0:sum(idx)-1)*spf/60;
            
            imagesc(tvecMin,zvecMic,AreaAF);
            colormap(parula);
            colorbar
            ylabel('Z-depth (microns)','FontSize',fs)
            xlabel('Time (min)','FontSize',fs)
            set(gca,'FontSize',fs)
            title(num2str(int),'FontSize',fs)
            maxL = max(abs(AreaAF),[],'all');
            minL = min(abs(AreaAF),[],'all');
            clim(gca,[minL, maxL]);
            pause;
            

            
        end
        



    end % T2 transition
end % movies

cd(od) % return to original directory
end % function