function [] = AreaRateHeatmap_LS_Rho_Wire(data,MovieNums,startCell)
% AreaRateHeatmap_LS  Produces a heatmap of the cell area rate
% This function was used to generate figure panels for the manuscript
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
% 04/12/23 Katie Bond



od = cd; % get original directory

fs = 16; % font size for figures

% length scale factor for X-Y (microns per pixel)
lsf = 0.163;

% length sacle factor for Z (mircons per pixel)
zsf = 1;
zvec = 4:19; % z-layers to select from z-stack
%zvecMic = zvec*zsf;
zvecMic = 1:16;

% deltaTSec is the time interval in seconds for caluclating the rate DA/DT
deltaTSec = 30;

% set the time window to choose from
tWindow = 1:250;


% use the function b2r to generate a colormap, flip so it goes from red to
% blue
b2rcmap = flip(b2r(-6,6),1);



Nmovies = numel(MovieNums);
for mn = 1:Nmovies
    
    % Move to directory and load required data
    cd(data(MovieNums(mn)).Source)
    %cd ..

    spf = data.SecPerFrame; % seconds per frame conversion factor
    
    % caluclate delta T for frames
    deltaTframes = round(deltaTSec/spf);


    % Load cell Area matrix and convert to microns^2
    loadGdata = load('GeometryData.mat');
    %AreaDataMat = loadGdata.Area(:,:,zvec)*lsf^2;
    AreaDataMat = loadGdata.Area(:,tWindow,zvec)*lsf^2;

    NumCells = size(AreaDataMat,1);

    for int=startCell:NumCells

        CellArea = squeeze(AreaDataMat(int,:,:))';
        finite = isfinite(CellArea);
        startFrame = find(all(finite,1),1,'first');
        stopFrame   = find(all(finite,1),1,'last');


        tvec = startFrame:stopFrame;
        
        if isempty(tvec)
            continue;
        end


        AreaA = CellArea(:,tvec);

      
        % smooth out the area matrix
        AreaAF = filterImage3DpaddedEdges(AreaA, 'Gauss', 2);

        % calculate the rate of change and convert to microns^2/min
        [dAdtmat] = RateOverDeltaT( AreaAF, deltaTframes);
        dAdtmat = dAdtmat/spf*60;

        tvecMin = (1:size(dAdtmat,2))*spf/60;
        
        if isempty(tvecMin) || tvecMin(end) < 5
            continue;
        end
        imagesc(tvecMin,zvecMic,dAdtmat);
        colormap(b2rcmap);
        caxis(gca,[-1 1]);
        colorbar
        ylabel('Z-depth (microns)','FontSize',fs)
        xlabel('Time (min)','FontSize',fs)
        set(gca,'FontSize',fs)
        title(num2str(int),'FontSize',fs)
        maxL = max(abs(dAdtmat),[],'all');
        caxis(gca,[-maxL, maxL]);
        
        
        choice = menu('Choose one: ','Select time window','Skip');
        if choice ==1
            StartStop = input('Input start and stop vec in min: ');
            
            
            idx = tvecMin>=StartStop(1) & tvecMin <= StartStop(2);
            
            
            dAdtmat = dAdtmat(:,idx);

            
            tvecMin = (0:sum(idx)-1)*spf/60;
            
            imagesc(tvecMin,zvecMic,dAdtmat);
            colormap(b2rcmap);
            caxis(gca,[-1 1]);
            colorbar
            ylabel('Z-depth (microns)','FontSize',fs)
            xlabel('Time (min)','FontSize',fs)
            set(gca,'FontSize',fs)
            title(num2str(int),'FontSize',fs)
            maxL = max(abs(dAdtmat),[],'all');
            %caxis(gca,[-maxL, maxL]);
            caxis(gca,[-10 10]);
            pause;
            

            
        end
        



    end % T2 transition
end % movies

cd(od) % return to original directory
end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                             SUBFUNCTIONS
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rate] = RateOverDeltaT( Signal2D, shift)
% INPUT:   
%           Signal2D: 2D matrix of cell area where rows are depth and
%           columns are time.
%           shift = Shift in number of frames to take the difference.


sh=shift;
tlen = size(Signal2D,2);

% cropped trajectory left
mat1 = Signal2D(:,1:(tlen-sh),:);
% cropped trajectory right
mat2 = Signal2D(:,(1+sh):tlen,:);
% gradient is difference between the two cropped trajectories divided by
% the time shift and has units of frames^-1.
Rate = (mat2-mat1)/sh;

end % of function
