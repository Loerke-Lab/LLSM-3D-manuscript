function [] = MyosinRateHeatmap_LS(data,MovieNums,startCell)
% MyosinRateHeatmap_LS  Produces a heatmap of the cell myosin rate
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
% 04/01/22 Tim Vanderleest



od = cd; % get original directory

fs = 16; % font size for figures

% length scale factor for X-Y (microns per pixel)
lsf = 0.104;

% length scale factor for Z (mircons per pixel)
zsf = 0.2635;
zvec = 1:60; % z-layers to select from z-stack
zvecMic = zvec*zsf;

% deltaTSec is the time interval in seconds for caluclating the rate DA/DT
deltaTSec = 30; 


% use the function b2r to generate a colormap, flip so it goes from red to
% blue
b2rcmap = flip(b2r(-6,6),1);
b2rcmap2 = flip(b2r2(-6,6),1);



Nmovies = numel(MovieNums);
for mn = 1:Nmovies
    
    % Move to directory and load required data
    cd(data(MovieNums(mn)).Source)
    cd ..

    spf = data.SecPerFrame; % seconds per frame conversion factor
    
    % caluclate delta T for frames
    deltaTframes = round(deltaTSec/spf);


    % % Load cell Area matrix and convert to microns^2
    % loadGdata = load('GeometryData.mat');
    % AreaDataMat = loadGdata.Area(:,:,zvec)*lsf^2;
    % 
    % NumCells = size(AreaDataMat,1);

    % load myosin intensity array (interior - M for medial) Myosin in channel 1
    loadMyo = load('CellBasedIntensityData.mat');
    myoArrayI = loadMyo.CellBasedIntensityData.InteriorIntensity_Myo(:,:,zvec);
    myoArrayB = loadMyo.CellBasedIntensityData.BoundaryIntensity_Myo(:,:,zvec);
    %myoArray = loadMyo.CellBasedIntensityData.JunctionIntensity_Myo(:,:,zvec);

    myoArray = myoArrayI + myoArrayB;

    NumCells = size(myoArray,1); % number of cells to loop over


    % Load cell Area matrix and convert to microns^2
    loadGdata = load('GeometryData.mat');
    AreaDataMat = loadGdata.Area(:,:,zvec)*lsf^2;

    for int=startCell:NumCells

        CellMyo = squeeze(myoArray(int,:,:))';
        finite = isfinite(CellMyo);
        startFrame = find(all(finite,1),1,'first');
        stopFrame   = find(all(finite,1),1,'last');


        tvec = startFrame:stopFrame;
        
        if isempty(tvec)
            continue;
        end


        MyoA = CellMyo(:,tvec);
      
        % smooth out the area matrix
        MyoAF = filterImage3DpaddedEdges(MyoA, 'Gauss', 2);

        % calculate the rate of change and convert to microns^2/min
        [dMdtmat] = RateOverDeltaT( MyoAF, deltaTframes);
        dMdtmat = dMdtmat/spf*60;


        tvecMin = (1:size(dMdtmat,2))*spf/60;
        
        if isempty(tvecMin) || tvecMin(end) < 5
            continue;
        end



        % add area data KL 3/19/24
        CellArea = squeeze(AreaDataMat(int,:,:))';
        AreaA = CellArea(:,tvec);

        % smooth out the area matrix
        AreaAF = filterImage3DpaddedEdges(AreaA, 'Gauss', 2);

        % calculate the rate of change and convert to microns^2/min
        [dAdtmat] = RateOverDeltaT( AreaAF, deltaTframes);
        dAdtmat = dAdtmat/spf*60;

        ax(1) = subplot(1,2,1);
        imagesc(tvecMin,zvecMic,dAdtmat);
        %colormap(b2rcmap2);
        caxis(gca,[-1 1]);
        colorbar
        ylabel('Z-depth (microns)','FontSize',fs)
        xlabel('Time (min)','FontSize',fs)
        set(gca,'FontSize',fs)
        title('Area Rate',num2str(int),'FontSize',fs)
        maxL = max(abs(dAdtmat),[],'all');
        caxis(gca,[-maxL, maxL]);

        ax(2) = subplot(1,2,2);
        imagesc(tvecMin,zvecMic,dMdtmat);
        %colormap(b2rcmap);
        caxis(gca,[-1 1]);
        colorbar
        ylabel('Z-depth (microns)','FontSize',fs)
        xlabel('Time (min)','FontSize',fs)
        set(gca,'FontSize',fs)
        title('Myosin Rate',num2str(int),'FontSize',fs)
        maxL = max(abs(dMdtmat),[],'all');
        caxis(gca,[-maxL, maxL]);

        colormap(ax(1),b2rcmap2);
        colormap(ax(2),b2rcmap);
        
        
        choice = menu('Choose one: ','Select time window','Skip');
        if choice ==1
            StartStop = input('Input start and stop vec in min: ');
            
            
            idx = tvecMin>=StartStop(1) & tvecMin <= StartStop(2);
            
            
            dMdtmat = dMdtmat(:,idx);
            dAdtmat = dAdtmat(:,idx);

            
            tvecMin = (0:sum(idx)-1)*spf/60;
            
            ax(1) = subplot(1,2,1);
            imagesc(tvecMin,zvecMic,dAdtmat);
            %colormap(b2rcmap2);
            caxis(gca,[-1 1]);
            colorbar
            ylabel('Z-depth (microns)','FontSize',fs)
            xlabel('Time (min)','FontSize',fs)
            set(gca,'FontSize',fs)
            title('Area Rate',num2str(int),'FontSize',fs)
            maxL = max(abs(dAdtmat),[],'all');
            caxis(gca,[-maxL, maxL]);
    
            ax(2) = subplot(1,2,2);
            imagesc(tvecMin,zvecMic,dMdtmat);
            %colormap(b2rcmap);
            caxis(gca,[-1 1]);
            colorbar
            ylabel('Z-depth (microns)','FontSize',fs)
            xlabel('Time (min)','FontSize',fs)
            set(gca,'FontSize',fs)
            title('Myosin Rate',num2str(int),'FontSize',fs)
            maxL = max(abs(dMdtmat),[],'all');
            caxis(gca,[-maxL, maxL]);
    
            colormap(ax(1),b2rcmap2);
            colormap(ax(2),b2rcmap);
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


function newmap = b2r(cmin_input,cmax_input)
%BLUEWHITERED   Blue, white, and red color map.
%   this matlab file is designed to draw anomaly figures, the color of
%   the colorbar is from blue to white and then to red, corresponding to 
%   the anomaly values from negative to zero to positive, respectively. 
%   The color white always correspondes to value zero. 
%   
%   You should input two values like caxis in matlab, that is the min and
%   the max value of color values designed.  e.g. colormap(b2r(-3,5))
%   
%   the brightness of blue and red will change according to your setting,
%   so that the brightness of the color corresponded to the color of his
%   opposite number
%   e.g. colormap(b2r(-3,6))   is from light blue to deep red
%   e.g. colormap(b2r(-3,3))   is from deep blue to deep red
%
%   I'd advise you to use colorbar first to make sure the caxis' cmax and cmin
%
%   by Cunjie Zhang
%   2011-3-14
%   find bugs ====> email : daisy19880411@126.com
%  
%   Examples:
%   ------------------------------
%   figure
%   peaks;
%   colormap(b2r(-6,8)), colorbar, title('b2r')
%   


%% check the input
if nargin ~= 2
   disp('input error');
   disp('input two variables, the range of caxis , for example : colormap(b2r(-3,3))')
end

if cmin_input >= cmax_input
    disp('input error')
    disp('the color range must be from a smaller one to a larger one')
end

%% control the figure caxis 
lims = get(gca, 'CLim');   % get figure caxis formation
caxis([cmin_input cmax_input])

%% color configuration : from blue to light blue to white untill to red

% red_top     = [1 0 0];
% white_middle= [1 1 1];
% blue_bottom = [0 0 1];
red_top     = [0.4940 0.1840 0.5560];
white_middle= [1 1 1];
blue_bottom = [0.4660 0.6740 0.1880];

%% color interpolation 

color_num = 250;   
color_input = [blue_bottom;  white_middle;  red_top];
oldsteps = linspace(-1, 1, length(color_input));
newsteps = linspace(-1, 1, color_num);  

%% Category Discussion according to the cmin and cmax input

%  the color data will be remaped to color range from -max(abs(cmin_input,abs(cmax_input)))
%  to max(abs(cmin_input,abs(cmax_input))) , and then squeeze the color
%  data in order to make suere the blue and red color selected corresponded
%  to their math values

%  for example :
%  if b2r(-3,6) ,the color range is from light blue to deep red ,
%  so that the color at -3 is light blue while the color at 3 is light red
%  corresponded

%% Category Discussion according to the cmin and cmax input
% first : from negative to positive
% then  : from positive to positive
% last  : from negative to negative

newmap_all = NaN(size(newsteps,2),3);
if (cmin_input < 0)  &&  (cmax_input > 0)
    
    
    if abs(cmin_input) < cmax_input 
         
        % |--------|---------|--------------------|    
      % -cmax      cmin       0                  cmax         [cmin,cmax]
      %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))
 
       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       start_point = round((cmin_input+cmax_input)/2/cmax_input*color_num);
       newmap = squeeze(newmap_all(start_point:color_num,:));
       
    elseif abs(cmin_input) >= cmax_input
        
         % |------------------|------|--------------|    
       %  cmin                0     cmax          -cmin         [cmin,cmax]
       %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))       
       
       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       end_point = round((cmax_input-cmin_input)/2/abs(cmin_input)*color_num);
       newmap = squeeze(newmap_all(1:end_point,:));
    end
    
       
elseif cmin_input >= 0

       if lims(1) < 0 
           disp('caution:')
           disp('there are still values smaller than 0, but cmin is larger than 0.')
           disp('some area will be in red color while it should be in blue color')
       end
        % |-----------------|-------|-------------|    
      % -cmax               0      cmin          cmax         [cmin,cmax]
      %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))
 
       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       start_point = round((cmin_input+cmax_input)/2/cmax_input*color_num);
       newmap = squeeze(newmap_all(start_point:color_num,:));

elseif cmax_input <= 0

       if lims(2) > 0 
           disp('caution:')
           disp('there are still values larger than 0, but cmax is smaller than 0.')
           disp('some area will be in blue color while it should be in red color')
       end
       
         % |------------|------|--------------------|    
       %  cmin         cmax    0                  -cmin         [cmin,cmax]
       %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))       

       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       end_point = round((cmax_input-cmin_input)/2/abs(cmin_input)*color_num);
       newmap = squeeze(newmap_all(1:end_point,:));
end
    

end


function newmap = b2r2(cmin_input,cmax_input)
%BLUEWHITERED   Blue, white, and red color map.
%   this matlab file is designed to draw anomaly figures, the color of
%   the colorbar is from blue to white and then to red, corresponding to 
%   the anomaly values from negative to zero to positive, respectively. 
%   The color white always correspondes to value zero. 
%   
%   You should input two values like caxis in matlab, that is the min and
%   the max value of color values designed.  e.g. colormap(b2r(-3,5))
%   
%   the brightness of blue and red will change according to your setting,
%   so that the brightness of the color corresponded to the color of his
%   opposite number
%   e.g. colormap(b2r(-3,6))   is from light blue to deep red
%   e.g. colormap(b2r(-3,3))   is from deep blue to deep red
%
%   I'd advise you to use colorbar first to make sure the caxis' cmax and cmin
%
%   by Cunjie Zhang
%   2011-3-14
%   find bugs ====> email : daisy19880411@126.com
%  
%   Examples:
%   ------------------------------
%   figure
%   peaks;
%   colormap(b2r(-6,8)), colorbar, title('b2r')
%   


%% check the input
if nargin ~= 2
   disp('input error');
   disp('input two variables, the range of caxis , for example : colormap(b2r(-3,3))')
end

if cmin_input >= cmax_input
    disp('input error')
    disp('the color range must be from a smaller one to a larger one')
end

%% control the figure caxis 
lims = get(gca, 'CLim');   % get figure caxis formation
caxis([cmin_input cmax_input])

%% color configuration : from blue to light blue to white untill to red

red_top     = [1 0 0];
white_middle= [1 1 1];
blue_bottom = [0 0 1];

%% color interpolation 

color_num = 250;   
color_input = [blue_bottom;  white_middle;  red_top];
oldsteps = linspace(-1, 1, length(color_input));
newsteps = linspace(-1, 1, color_num);  

%% Category Discussion according to the cmin and cmax input

%  the color data will be remaped to color range from -max(abs(cmin_input,abs(cmax_input)))
%  to max(abs(cmin_input,abs(cmax_input))) , and then squeeze the color
%  data in order to make suere the blue and red color selected corresponded
%  to their math values

%  for example :
%  if b2r(-3,6) ,the color range is from light blue to deep red ,
%  so that the color at -3 is light blue while the color at 3 is light red
%  corresponded

%% Category Discussion according to the cmin and cmax input
% first : from negative to positive
% then  : from positive to positive
% last  : from negative to negative

newmap_all = NaN(size(newsteps,2),3);
if (cmin_input < 0)  &&  (cmax_input > 0)
    
    
    if abs(cmin_input) < cmax_input 
         
        % |--------|---------|--------------------|    
      % -cmax      cmin       0                  cmax         [cmin,cmax]
      %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))
 
       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       start_point = round((cmin_input+cmax_input)/2/cmax_input*color_num);
       newmap = squeeze(newmap_all(start_point:color_num,:));
       
    elseif abs(cmin_input) >= cmax_input
        
         % |------------------|------|--------------|    
       %  cmin                0     cmax          -cmin         [cmin,cmax]
       %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))       
       
       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       end_point = round((cmax_input-cmin_input)/2/abs(cmin_input)*color_num);
       newmap = squeeze(newmap_all(1:end_point,:));
    end
    
       
elseif cmin_input >= 0

       if lims(1) < 0 
           disp('caution:')
           disp('there are still values smaller than 0, but cmin is larger than 0.')
           disp('some area will be in red color while it should be in blue color')
       end
        % |-----------------|-------|-------------|    
      % -cmax               0      cmin          cmax         [cmin,cmax]
      %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))
 
       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       start_point = round((cmin_input+cmax_input)/2/cmax_input*color_num);
       newmap = squeeze(newmap_all(start_point:color_num,:));

elseif cmax_input <= 0

       if lims(2) > 0 
           disp('caution:')
           disp('there are still values larger than 0, but cmax is smaller than 0.')
           disp('some area will be in blue color while it should be in red color')
       end
       
         % |------------|------|--------------------|    
       %  cmin         cmax    0                  -cmin         [cmin,cmax]
       %    squeeze(colormap(round((cmin+cmax)/2/cmax),size(colormap)))       

       for j=1:3
           newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
       end
       end_point = round((cmax_input-cmin_input)/2/abs(cmin_input)*color_num);
       newmap = squeeze(newmap_all(1:end_point,:));
end
    

end

