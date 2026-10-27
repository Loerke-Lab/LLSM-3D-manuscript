function [myosinAutoCorr,maxLag] = myosinZLayers_LS(data,MovieNum,direction)

% this function runs cross correlations of myosin from the apical layer 
% of the germband down consecutive layers to the most basal layer.

% INPUT:    data: a data structure containing information about the movie
%                 including Source, ImageFileList, SecPerFrame
%           MovieNum: the movie number to be analyzed if there are several
%                     movies stored in the data structure (will be 1 if there is only
%                     1 movie)
%           direction: the direction to run the cross correlation (i.e.
%           'AtoB', 'BtoA', 'MtoB', 'MtoA')

% OUTPUT:   lagsCellAndZ: a 3D matrix consisting of lags for cross correlations
%                         for cells at all time points and z layers
%                         rows: individual cells
%                         columns: cross correlation lags in minutes
%                         3rd dimension: z layers
%           XCorrCellAndZ: a 3D matrix consisting of cross correlation
%                          results for all cells at all time points and z layers
%                          rows: individual cells
%                          columns: normalized cross correlation values
%                          3rd dimension: z layers


% get original directory
od = cd;

% set sigma for interpolation
sigma = 1;

% set lag for cross correlation
maxLag = 20;

% set depth spacing for light sheet data
Zsf = 0.263;

zvec = 1:60; % z-layers to process
zvecMic = zvec*Zsf;

% initialize arrays to save all movies data
NlayersAll = [];
avgCCNormAll = {};
avgXCAllMovies = [];
myosinAutoCorr = [];
countAll = 1;

shiftSec = 60;  % delta T for computing the rate
spf = data(1).SecPerFrame;
shiftF = round(shiftSec/spf);

% loop over the movies
Nmovies = numel(MovieNum);

%figure(1)
for mn = 1:Nmovies

    % set the movie number to use
    mnToUse = MovieNum(mn);

    % get seconds per frame for the movie
    spf = data(mnToUse).SecPerFrame;
    
    % move to directory with the movie data
    cd(data(mnToUse).Source)
    cd ..

    % get intensity data
    loadMyo = load('CellBasedIntensityData.mat').CellBasedIntensityData;
    %myosinArray = loadMyo.BoundaryIntensity_Myo(:,:,zvec);
    myosinArray = loadMyo.InteriorIntensity_Myo(:,:,zvec);


    % get number of cells in the array
    Ncells = size(myosinArray,1);

    % get number of layers in the array
    Nlayers = size(myosinArray,3);
    
    % initialize results matrices
    lagsCellAndZ = zeros(Ncells,(maxLag*2)+1,Nlayers);
    XCorrCellAndZ = zeros(Ncells,(maxLag*2)+1,Nlayers);


    count = 1;

    % get the middle layer for propagations
    middleLayer = round(Nlayers/2);


    % set variables for directionality of cross correlation
    switch direction
        case 'AtoB'
            forLoop = 1:1:Nlayers;
            layerStart = 1;
            cmap = turbo(Nlayers);
            maxDepth = (Nlayers*Zsf)-Zsf;
            depthMat = [0:Zsf:maxDepth];
            txtHM = ['Cross Correlation Apical to Basal',' (Movie ',num2str(MovieNum),')'];
        case 'BtoA'
            forLoop = Nlayers:-1:1;
            layerStart = Nlayers;
            cmap = turbo(Nlayers);
            maxDepth = (Nlayers*Zsf)-Zsf;
            depthMat = [maxDepth:-Zsf:0];
            txtHM = ['Cross Correlation Basal to Apical',' (Movie ',num2str(MovieNum),')'];
        case 'MtoA'
            forLoop = middleLayer:-1:1;
            layerStart = middleLayer;
            cmap = turbo(middleLayer+1);
            middleDepth = (middleLayer*Zsf)-Zsf;
            depthMat = [middleDepth:-Zsf:0];
            txtHM = ['Cross Correlation Middle to Apical',' (Movie ',num2str(MovieNum),')'];
        case 'MtoB'
            forLoop = middleLayer:1:Nlayers;
            layerStart = middleLayer;
            cmap = turbo(middleLayer+1);
            middleDepth = (middleLayer*Zsf)-Zsf;
            maxDepth = (Nlayers*Zsf)-Zsf;
            depthMat = [middleDepth:Zsf:maxDepth];
            txtHM = ['Cross Correlation Middle to Basal',' (Movie ',num2str(MovieNum),')'];
    end
    


    % loop over the layers
    for j = forLoop


         % set the z layer to compare to (inside the loop)
        layerCompare = j;

        % loop over the number of cells in the movie
        for i = 1:Ncells
    
            cellMyosinArrayNorm = [];

            % get the cell ID (inside the loop)
            %cellID = cellNumbers(i,:);

            % find a time vector in which the junction is present at all
            % depths, if it is empty of has a duration of less than 3 minutes
            % skip to next junction
            cellMyosinArray = squeeze(myosinArray(i,:,:))'; % rows are layers, columns are time points

            [~, tvecF] = FiniteSequenceFromNans(mean(cellMyosinArray,1));%,1);
            if isempty(tvecF) || length(tvecF)*spf/60 <3
                continue;
            end

            cellMyosinArray = cellMyosinArray(:,tvecF);

            % compute rates of cell area change
            [cellMyosinArray] = ShiftedRate( cellMyosinArray, shiftF);

            % get the start myosin layer
            myosinLayerStart = cellMyosinArray(layerStart,:);
    
            % get the myosin layer to compare to
            myosinLayerCompare = cellMyosinArray(layerCompare,:);

            % normalize the entire matrix (TV method)
            normMethod = 1;
            Nrows = size(cellMyosinArray,1);    
    
            if normMethod == 1
                % normalize the entire matrix, 
                cellMyosinArrayNorm = (cellMyosinArray - mean(cellMyosinArray,'all'))/std(cellMyosinArray,[],'all');

                % get the start myosin layer
                myosinLayerStartNorm = cellMyosinArrayNorm(layerStart,:);
    
                % get the myosin layer to compare to
                myosinLayerCompareNorm = cellMyosinArrayNorm(layerCompare,:);
            end
            

            if normMethod == 2
                myosinLayerStartNorm = (myosinLayerStart - mean(myosinLayerStart,'all'))/std(myosinLayerStart,[],'all');
                myosinLayerCompareNorm = (myosinLayerCompare - mean(myosinLayerCompare,'all'))/std(myosinLayerCompare,[],'all');

                for l = forLoop
                    cellMyosinArrayNorm(l,:) = (cellMyosinArray(l,:) - mean(cellMyosinArray(l,:),'all'))/std(cellMyosinArray(l,:),[],'all');
                end
            end

            if normMethod == 3
            end

            if normMethod == 4

            numFrames = size(myosinLayerStart,2);
            frameVec = (1:numFrames);
            timeVec = (frameVec*spf)/60;

            [myosinLayerStartNorm] = interpolateAreaOrMyosin(myosinLayerStart,sigma,timeVec,numFrames);
            [myosinLayerCompareNorm] = interpolateAreaOrMyosin(myosinLayerCompare,sigma,timeVec,numFrames);

                % for l = forLoop
                %     cellMyosinArrayNorm(l,:) = interpolateAreaOrMyosin(cellMyosinArray(l,:),sigma,timeVec,numFrames);
                % end

            end


            % % get the time vector for plotting and interpolation (area start layer)
            % numFrames = size(myosinLayerStart,2);
            % frameVec = (1:numFrames);
            % timeVec = (frameVec*spf)/60;
            % 
            % % skip over the cell if the time vector is empty
            % if isempty(timeVec)
            %     continue
            % end
            % 
            % % run the interpolation for the most apical layer
            % [myosinLayerStartInt] = interpolateAreaOrMyosin(myosinLayerStart,sigma,timeVec,numFrames);
            % 
            % % get the time vector for plotting and interpolation (area layer to compare)
            % numFrames = size(myosinLayerCompare,2);
            % frameVec = (1:numFrames);
            % timeVec = (frameVec*spf)/60;
            % 
            % % skip over the cell if the time vector is empty
            % if isempty(timeVec)
            %     continue
            % end
            % 
            % % run the interpolation for the layer to compare
            % [myosinLayerCompareInt] = interpolateAreaOrMyosin(myosinLayerCompare,sigma,timeVec,numFrames);
    

            % cross correlate the area layers
            %[c,lags] = xcorr(myosinLayerStart,myosinLayerCompare,maxLag,'normalized');
            [c,lags] = xcorr(myosinLayerStartNorm,myosinLayerCompareNorm,maxLag,'normalized');
            lagsTime = (lags*spf)/60;

            % save cross correlation results for each cell at each layer
            lagsCellAndZ(i,:,j) = lagsTime;
            XCorrCellAndZ(i,:,j) = c;

            % subplot(2,2,1)
            % imagesc(tvecF*spf/60,zvecMic,cellMyosinArray);colorbar
            % title('Myosin')
            % 
            % subplot(2,2,3)
            % imagesc(tvecF*spf/60,zvecMic,cellMyosinArrayNorm);colorbar
            % title('Myosin Norm')

        end

        % find out which rows are all zeros.
        ZeroRows = all(XCorrCellAndZ(:,:,j) == 0, 2);
        % assign all columns of the all-zero rows to nans.
        XCorrCellAndZ(ZeroRows,:,j) = NaN;
        lagsCellAndZ(ZeroRows,:,j) = NaN;
        
        % average cross correlation
        avgXC(j,:) = mean(XCorrCellAndZ(:,:,j),1,'omitnan');
        avgLags(j,:) = mean(lagsCellAndZ(:,:,j),1,'omitnan');


        % % plot the average cross correlation
        % plot(avgLags(j,:),avgXC(j,:),'LineWidth',2,'Color',cmap(count,:));
        % set(gca,'FontSize',16)
        % txt = ['Movie ',num2str(MovieNum),' (blue to red in layers)'];%' Z depth ',num2str(zSpacing*(j-1)),' microns'];
        % title(txt)
        % xlabel('Lag (minutes)')
        % ylabel('Correlation')
        % grid on
        % %ylim([-40 150])
        % %pause
        % hold on
        % set(gca,'Fontsize',20)
        % 
        count = count + 1;

    end
    
   %  yLayers = repmat([1:Nlayers]',1,size(avgLags,2));
   % 
   %  depthMat = round(depthMat,2);
   %  depthLayers = repmat(depthMat',1,size(avgLags,2));
   % 
   %  % Convert each number in the array into a string
   %  CustomYLabels = string(depthMat);
   %  % depth divide
   %  dDivide = zSpacing*5;
   %  % Replace all but the fifth elements by spaces
   %  CustomYLabels(mod(depthMat,dDivide) ~= 0) = " ";
   % 
   %  avgXCMiddle = avgXC(middleLayer:end,:);
   % 
   %  % set the average array depending on direction
   %  switch direction
   %      case 'BtoA'
   %          avgXC = flipud(avgXC);
   %      case 'MtoA'
   %          avgXC = flipud(avgXC);
   %      case 'MtoB'
   %          avgXC = avgXCMiddle;
   %  end
   % 
   % 
   %  % plot a heatmap of the average cross correlation for all cells
   %  figure
   %  heatmap(avgXC,'Colormap',turbo,'ColorLimits',[-0.4 1]);%b2r(-0.25,1))
   %  %heatmap(flipud(avgXCMiddle),'Colormap',turbo,'ColorLimits',[-0.4 1]);%b2r(-0.25,1))
   %  Ax = gca;
   %  Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
   % % Ax.XDisplayLabels = CustomXLabels;
   %  Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
   %  Ax.YDisplayLabels = CustomYLabels;
   %  %txtHM = ['Cross Correlation Apical to Basal',' (Movie ',num2str(MovieNum),')'];
   %  set(gca,'FontSize',14)
   %  %sgtitle(txt,'FontSize',18) 
   %  title(txtHM)
   %  xlabel('Lag')
   %  ylabel('Depth (microns)')

    myosinAutoCorr = cat(1,myosinAutoCorr,XCorrCellAndZ);
    avgXCAllMovies = cat(3,avgXCAllMovies,avgXC);

    % reinitialize arrays for averaging
    avgXC = [];
    avgXCNorm = [];
    avgLags = [];
    count = 1;


end

    lagvecSec = lags*spf;

    % subplot(2,2,[2 4])
    % imagesc(lagvecSec,zvecMic,squeeze(myosinAutoCorr(15,:,:))');colorbar
    % title('Acorr')


    avgXCAllMovies = mean(avgXCAllMovies,3);

    figure
    for l = forLoop

        % plot the average cross correlation
        plot(lags*spf,avgXCAllMovies(l,:),'LineWidth',2,'Color',cmap(countAll,:));
        set(gca,'FontSize',16)
        txt = ['Myosin Rate AutoCorr (blue to red in layers)'];%' Z depth ',num2str(zSpacing*(j-1)),' microns'];
        title(txt)
        xlabel('Lag (seconds)')
        ylabel('Correlation')
        grid on
        %ylim([-40 150])
        %pause
        hold on
        set(gca,'Fontsize',20)
    
        countAll = countAll + 1;

    end

    lagvecSec = lags*spf;

    figure
    %subplot(1,2,1)
    imagesc(lagvecSec,zvecMic,avgXCAllMovies)
    xlabel('Lag (sec)')
    ylabel('Depth (microns)')
    title('Myosin Rate Autocorrelation')
    colorbar
    caxis([-0.4 1])
    set(gca,'FontSize',20)

    % subplot(1,2,2)
    % apical_zvec = round(5/Zsf);
    % lateral_zvec = round(10/Zsf);
    % basal_zvec  = round(15/Zsf);
    % plot(lagvecSec,mean(avgXCAllMovies(apical_zvec,:),1),'LineWidth',2);hold on
    % plot(lagvecSec,mean(avgXCAllMovies(lateral_zvec,:),1),'LineWidth',2);hold on
    % plot(lagvecSec,mean(avgXCAllMovies(basal_zvec,:),1),'LineWidth',2);hold on
    % xlabel('Lag (sec)')
    % ylabel('Cross-Correlation coeff.')
    % legend('5 microns','10 microns','15 microns')
    % set(gca,'FontSize',20)
    % grid on
    % %set(gca,'Ylim',[-0.8 0.2])


cd(od)


end



%%% SUBFUNCTIONS %%%

function [VectorOut, idx_out] = FiniteSequenceFromNans(VectorIn)
% This function takes a vector, which may contain NaN values throughout,
% and outputs the longest continuous sequence which contains no NaNs. It 
% also outputs the indices of that sequence. For example, the input 
% [ NaN NaN 1 4 NaN 7 9 12 NaN] gives the out put [7 9 12] with indices
% [6 7 8].

if all(isnan(VectorIn(:)))  % These 5 lines were added to return empty outputs
    VectorOut = [];         % when the input vector is all nans
    idx_out = [];
    return;
end


finiteIdx = find(isfinite(VectorIn));
x = [0 cumsum(diff(finiteIdx)~=1)];
idx_out = finiteIdx(x==mode(x));
VectorOut = VectorIn(idx_out);
    
end % end FiniteSequenceFromNans

function [outputArray] = interpolateAreaOrMyosin(inputArray,sigma,timeVec,numFrames)

% set the time bin for interpolation
maxTime = max(timeVec);
minTime = min(timeVec);

bins = linspace(minTime,maxTime,numFrames);

% set the value matrix and run the interpolation
valueMatrix = [timeVec(:),inputArray(:)];
[interpArrayVTime] = interpolatePointVals2Vec_gauss_wstd(valueMatrix(:,1:2), bins, sigma, 0 );

% pull out interpolated area values
arrayInt = interpArrayVTime(:,2);

% rearrange AreaArray back to original vector shape
ArrayVecInt = arrayInt';

% subtract interpolation vector from original vector
ArrayVecDetrend = inputArray - ArrayVecInt;

% normalize by subtracting the mean and dividing by the standard deviation
outputArray = (ArrayVecDetrend-mean(ArrayVecDetrend))/(std(ArrayVecDetrend));

end


function [Signalsh] = ShiftedRate( Signal, shift)
% INPUT:    vert1 = vertex trajectory 1
%           vert2 = vertex trajectory 2
%           shiftvec = vector with desired shifts for gradient (in frames)


sh=shift;
tlen = size(Signal,2);

% cropped trajectory left
mat1 = Signal(:,1:(tlen-sh));
% cropped trajectory right
mat2 = Signal(:,(1+sh):tlen);
% gradient is difference between the two cropped trajectories divided by
% the time shift and has units of frames^-1.
Signalsh = (mat2-mat1)/sh;

end % of function









