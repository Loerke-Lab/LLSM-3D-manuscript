function [] = cellAreaMyosinXC_LS(data,MovieNum,showSamples)

% this function will compare cell area oscillations with myosin intensity at various depths in the 
% light sheet data

% INPUT:    data: a data structure containing information about the movie
%                 including Source, ImageFileList, SecPerFrame
%           MovieNum: the movie number to be analyzed if there are several
%                     movies stored in the data structure (will be 1 if there is only
%                     1 movie)
%
% OUTPUT:   no specific output, just plots the area, myosin, autocorrelation, and cross correlation 
%           for each cell. plots average cross correlation for each movie.


% record original directory (to return to it at the end)
od = cd;

% set sigmas for interpolation
sigma = 1;

% set max layers and step space through layers
zvec = 1:60;
step = 1;
%layersToAnalyze = [1 19 32 49 65];
%layersToAnalyze = [1 32 65];
Alayer = 1;
Llayer = 30;%32;
Blayer = 60;%65;

Zsf = 0.2635;
zvecMic = zvec*Zsf;

% set lag for cross correlation
maxLag = 20;
xcorr_resultsMat  = [];
%zeroLagXC = {};
Nlags = (maxLag*2)+1;
avgXCAllMovies = [];

shiftSec = 60;  % delta T for computing the rate
spf = data(1).SecPerFrame;
shiftF = round(shiftSec/spf);


% loop through the movies
Nmovies = numel(MovieNum);
for mn=1:Nmovies

    % get seconds per frame for the movie (& start/end GBE)
    spf = data(MovieNum(mn)).SecPerFrame;

    % move to the parent directory of the source of movie
    cd(data(MovieNum(mn)).Source)
    cd ..

    % load cell area array
    loadArea = load('GeometryData.mat');
    areaArray = loadArea.Area(:,:,zvec);
    cellNumbers = loadArea.CellNumbers;

    % load myosin intensity array (interior - M for medial) Myosin in channel 1
    loadMyo = load('CellBasedIntensityData.mat');
    myoArray = loadMyo.CellBasedIntensityData.InteriorIntensity_Myo(:,:,zvec);
    %myoArray = data(MovieNum(mn)).CellIntensityNorm(:,:,zvec);

    % get the number of cells in the array
    Ncells = size(areaArray,1);

    % get number of layers in the array
    Nlayers = size(areaArray,3);

    % initialize results matrices
    % lagsCellAndZ = zeros(Ncells,(maxLag*2)+1,Nlayers);
    % XCorrCellAndZ = zeros(Ncells,(maxLag*2)+1,Nlayers);

    lagsCellAndZ = NaN(Ncells,(maxLag*2)+1,Nlayers);
    XCorrCellAndZ = NaN(Ncells,(maxLag*2)+1,Nlayers);

    count = 1;
    cmap = turbo(Nlayers);

    % loop over the layers
    %for l = 1:step:Nlayers%layersToAnalyze
    %for l = layersToAnalyze

        % loop over the cells in the movie
        for c = 1:Ncells

            % get the cell ID (inside the loop)
            cellID = cellNumbers(c,:);

            % find a time vector in which the cell is present at all
            % depths, if it is empty of has a duration of less than 3 minutes
            % skip to next cell
            % cellAreaArray = squeeze(areaArray(c,:,1:Nlayers))'; % rows are layers, columns are time points
            % cellMyoArray = squeeze(myoArray(c,:,1:Nlayers))';    

            cellAreaArray = squeeze(areaArray(c,:,:))'; % rows are layers, columns are time points
            cellMyoArray = squeeze(myoArray(c,:,:))';  

            [~, tvecF] = FiniteSequenceFromNans(mean(cellAreaArray,1));%,1);
            if isempty(tvecF) || length(tvecF)*spf/60 <3
                continue;
            end

            % % get the cell and myosin arrays for this layer and the time vector
            % cellAreaArray = cellAreaArray(l,tvecF);
            % cellMyoArray = cellMyoArray(l,tvecF);
            % 
            % % get the time vector for plotting and interpolation (area start layer)
            % numFrames = size(cellAreaArray,2);
            % frameVec = (1:numFrames);
            % timeVec = (frameVec*spf)/60;
            % 
            % % skip over the cell if the time vector is empty
            % if isempty(timeVec)
            %     continue
            % end
            % 
            % % run the interpolation for the most apical layer
            % [cellAreaInt] = interpolateAreaOrMyosin(cellAreaArray,sigma,timeVec,numFrames);
            % 
            % % get the time vector for plotting and interpolation (area layer to compare)
            % numFrames = size(cellMyoArray,2);
            % frameVec = (1:numFrames);
            % timeVec = (frameVec*spf)/60;
            % 
            % % skip over the cell if the time vector is empty
            % if isempty(timeVec)
            %     continue
            % end
            % 
            % % run the interpolation for the layer to compare
            % [cellMyoInt] = interpolateAreaOrMyosin(cellMyoArray,sigma,timeVec,numFrames);
            % 
            % % cross correlate area against myosin
            % [xc,lags] = xcorr(cellAreaInt,cellMyoInt,maxLag,'normalized');
            % lagsTime = (lags*spf)/60;
            % 
            % % save cross correlation results for each cell at each layer
            % lagsCellAndZ(c,:,l) = lagsTime;
            % XCorrCellAndZ(c,:,l) = xc;


            % get the cell and myosin arrays for all layers and the time vector
            cellAreaArray = cellAreaArray(:,tvecF);
            cellMyoArray = cellMyoArray(:,tvecF);

            % compute rates of cell area change and medial myosin intensity
            [cellAreaRate] = ShiftedRate( cellAreaArray, shiftF);
            [cellMyoRate] = ShiftedRate( cellMyoArray, shiftF);

            % perform the cross-correlation of the rates and store into a layer
            % of a results matrix.
            % normMethod 1 is the whole matrix is normalized, normMethod 2 is each
            % z-layer is normalized
            normMethod = 1;
            %[xcorr_result,lags,areaNorm,myoNorm] = XcorrFunction(cellAreaArray,cellMyoArray,maxLag,normMethod,spf);
            [xcorr_result,lags,areaNorm,myoNorm] = XcorrFunction(cellAreaRate,cellMyoRate,maxLag,normMethod,spf);
            xcorr_resultsMat = cat(3,xcorr_resultsMat,xcorr_result);

            lagvec = -maxLag:maxLag;
            lagvecSec = lagvec*spf;

            if showSamples
            
                subplot(2,3,1)
                imagesc(tvecF*spf/60,zvecMic,cellAreaArray);colorbar
                title('Area')
                
                subplot(2,3,2)
                imagesc(tvecF*spf/60,zvecMic,cellMyoArray);colorbar
                title('Myosin')
                
                subplot(1,3,3)
                imagesc(lagvecSec,zvecMic,xcorr_result);colorbar
                title('Xcorr')

                subplot(2,3,4)
                imagesc(tvecF*spf/60,zvecMic,areaNorm);colorbar
                title('Area Norm')
                
                subplot(2,3,5)
                imagesc(tvecF*spf/60,zvecMic,myoNorm);colorbar
                title('Myosin Norm')
                
                pause
            
            end


            % [acArea,lags] = xcorr(cellAreaInt,maxLag,'normalized');
            % [acMyo,lags] = xcorr(cellMyoInt,maxLag,'normalized');


            % figure(1)
            % subplot(3,1,1)
            % plot(timeVec,cellAreaInt,'LineWidth',2)
            % hold on
            % plot(timeVec,cellMyoInt,'LineWidth',2)
            % title('Cell Area & Myosin (detrended)')
            % legend('Area','Myosin','Location','northoutside','Orientation','horizontal')
            % xlabel('Time')
            % ylabel('Amplitude')
            % grid on
            % set(gca,'FontSize',18)
            % hold off
            % 
            % subplot(3,1,2)
            % plot(lagsTime,xc,'LineWidth',2)
            % xlabel('Lag (minutes)')
            % ylabel('Correlation')
            % title('Area vs Myosin (Interior) Correlation')
            % grid on
            % set(gca,'FontSize',18)
            % 
            % subplot(3,1,3)
            % plot(lagsTime,acArea,'LineWidth',2)
            % hold on 
            % plot(lagsTime,acMyo,'LineWidth',2)
            % title('Autocorrelations')
            % legend('Area','Myosin','Location','northoutside','Orientation','horizontal')
            % xlabel('Lag (minutes)')
            % ylabel('Correlation')
            % grid on
            % set(gca,'FontSize',18)
            % hold off
            % % txt = ['Cell ',num2str(cellID)];
            % % sgtitle(txt,'FontSize',18)
            % pause

        end

        % % find out which rows are all zeros.
        % ZeroRows = all(XCorrCellAndZ(:,:,l) == 0, 2);
        % % assign all columns of the all-zero rows to nans.
        % XCorrCellAndZ(ZeroRows,:,l) = NaN;
        % lagsCellAndZ(ZeroRows,:,l) = NaN;

        % average cross correlation
        % avgXC(l,:) = mean(XCorrCellAndZ(:,:,l),1,'omitnan');
        % avgLags(l,:) = mean(lagsCellAndZ(:,:,l),1,'omitnan');


        %avgXC = mean(xcorr_resultsMat,3,'omitnan');
        avgXC  = NaN(Nlayers,Nlags);
        for l = 1:step:Nlayers

        XMLz = squeeze(xcorr_resultsMat(l,:,:))';
            
        avgXC(l,:) = mean(XMLz,'omitnan');

        % plot the average cross correlation
        plot(lags*spf/60,avgXC(l,:),'LineWidth',2,'Color',cmap(l,:));
        set(gca,'FontSize',16)
        txt = ['Movie ',num2str(MovieNum),' (blue to red in layers)'];%' Z depth ',num2str(zSpacing*(j-1)),' microns'];
        title(txt)
        xlabel('Lag (minutes)')
        ylabel('Correlation')
        grid on
        %ylim([-40 150])
        %pause
        hold on
        set(gca,'Fontsize',20)
        %legend('0 microns','4 microns','8 microns','12 microns','16 microns')
    
        count = count + step;
    
    
        end

    % make boxplot of xc value at zero lag
    apicalXC = XCorrCellAndZ(:,maxLag+1,Alayer);
    lateralXC = XCorrCellAndZ(:,maxLag+1,Llayer);
    basalXC = XCorrCellAndZ(:,maxLag+1,Blayer);

    apicalXC = apicalXC(~isnan(apicalXC));
    lateralXC = lateralXC(~isnan(lateralXC));
    basalXC = basalXC(~isnan(basalXC));

    y = [apicalXC; lateralXC; basalXC];
    numCells = size(apicalXC,1);
    g1 = repmat({'Apical'},numCells,1);
    g2 = repmat({'Lateral'},numCells,1);
    g3 = repmat({'Basal'},numCells,1);
    g = [g1; g2; g3];
    figure
    boxplot(y,g,'Symbol','')
    ylabel('XCorr at Zero Lag Frame')
    title(['Ncells=',num2str(numCells)])
    grid on
    %set(gca,'Ylim',[0 0.03])
    set(findobj(gca,'type','line'),'linew',2)
    set(gca,'FontSize',20)
    ylim([-0.9 0.7])

    y2 = [mean(apicalXC); mean(lateralXC); mean(basalXC)];
    hold on
    plot(y2,'dg','MarkerFaceColor','g','MarkerSize',10)

    
   %  % plot heatmap of cross correlation results
   %  maxDepth = (Nlayers*zSpacing)-zSpacing;
   %  depthMat = [0:zSpacing:maxDepth];
   %  depthMat = round(depthMat,2);
   % 
   % 
   %  % Convert each number in the array into a string
   %  CustomYLabels = string(depthMat);
   %  % depth divide
   %  dDivide = zSpacing*5;
   %  % Replace all but the fifth elements by spaces
   %  CustomYLabels(mod(depthMat,dDivide) ~= 0) = " ";
   % 
   % 
   %  % plot a heatmap of the average cross correlation for all cells
   %  figure
   %  heatmap(avgXC,'Colormap',turbo);%b2r(-0.25,1))
   %  Ax = gca;
   %  Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
   % % Ax.XDisplayLabels = CustomXLabels;
   %  Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
   %  Ax.YDisplayLabels = CustomYLabels;
   %  txtHM = ['Cross Correlation at Each Layer',' (Movie ',num2str(MovieNum),')'];
   %  set(gca,'FontSize',14)
   %  %sgtitle(txt,'FontSize',18) 
   %  title(txtHM)
   %  xlabel('Lag')
   %  ylabel('Depth (microns)')
   %  caxis([-0.6 0.2])
   %  set(gca,'Fontsize',20)


    avgXCAllMovies = cat(3,avgXCAllMovies,avgXC);

    zeroLagMovie = squeeze(xcorr_resultsMat(:,maxLag+1,:));
    zeroLagXC{:,mn} = zeroLagMovie;

    % reinitialize arrays for averaging
    avgXC = [];
    avgLags = [];
    count = 1;
    xcorr_resultsMat = [];

end

    avgXCAllMovies = mean(avgXCAllMovies,3,'omitnan');
    
   %  % plot heatmap of cross correlation results (for all movies)
   %  maxDepth = (Nlayers*Zsf)-Zsf;
   %  depthMat = [0:Zsf:maxDepth];
   %  depthMat = round(depthMat,2);
   % 
   % 
   %  % Convert each number in the array into a string
   %  CustomYLabels = string(depthMat);
   %  % depth divide
   %  dDivide = Zsf*5;
   %  % Replace all but the fifth elements by spaces
   %  CustomYLabels(mod(depthMat,dDivide) ~= 0) = " ";
   % 
   % 
   %  % plot a heatmap of the average cross correlation for all cells
   %  figure
   %  heatmap(avgXCAllMovies,'Colormap',turbo);%b2r(-0.25,1))
   %  Ax = gca;
   %  Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
   % % Ax.XDisplayLabels = CustomXLabels;
   %  Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
   %  Ax.YDisplayLabels = CustomYLabels;
   %  txtHM = ['Cross Correlation at Each Layer',' (Movie ',num2str(MovieNum),')'];
   %  set(gca,'FontSize',14)
   %  %sgtitle(txt,'FontSize',18) 
   %  title(txtHM)
   %  xlabel('Lag')
   %  ylabel('Depth (microns)')
   %  caxis([-0.6 0.2])
   %  set(gca,'Fontsize',20)


    lagvecSec = lags*spf;

    figure
    subplot(1,2,1)
    imagesc(lagvecSec,zvecMic,avgXCAllMovies)
    xlabel('Lag (sec)')
    ylabel('Depth (microns)')
    %title('Cross Sectional Area vs Medial Myosin')
    title('Area vs Medial Myosin Rates')
    colorbar
    caxis([-0.4 0.1])
    set(gca,'FontSize',20)

    subplot(1,2,2)
    %top_zvec = 1;
    apical_zvec = round(5/Zsf);
    lateral_zvec = round(10/Zsf);
    basal_zvec  = round(15/Zsf);


    % apical_zvec = 1:20;
    % lateral_zvec = 21:40;
    % basal_zvec  = 41:60;
    %plot(lagvecSec,mean(avgXCAllMovies(top_zvec,:),1),'LineWidth',2);hold on
    plot(lagvecSec,mean(avgXCAllMovies(apical_zvec,:),1),'LineWidth',2);hold on
    plot(lagvecSec,mean(avgXCAllMovies(lateral_zvec,:),1),'LineWidth',2);hold on
    plot(lagvecSec,mean(avgXCAllMovies(basal_zvec,:),1),'LineWidth',2);hold on
    xlabel('Lag (sec)')
    ylabel('Cross-Correlation coeff.')
    legend('5 microns','10 microns','15 microns')
    %legend('Apical','Lateral','Basal')
    set(gca,'FontSize',20)
    grid on
    set(gca,'Ylim',[-0.4 0.1])



    % get the average,standard deviation, and standard error
    avgXC = cellfun(@(zeroLagXC) mean(zeroLagXC,2), zeroLagXC, 'UniformOutput', false);
    avgXC = mean(cell2mat(avgXC),2);
    stdXC = cellfun(@(zeroLagXC) std(zeroLagXC,[],2), zeroLagXC, 'UniformOutput', false);
    stdXC = std(cell2mat(stdXC),[],2);
    nCells = cellfun(@size,zeroLagXC,'uni',false);
    nCells = cell2mat(nCells);
    nCells = sum(nCells(2:2:end));
    seXC = (stdXC/sqrt(nCells));


    % plot the values at zero lag with standard error envelope
    figure
    x = zvecMic';
    y = avgXC;
    e = seXC;
    lo = y - e;
    hi = y + e;
    hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], 'k');
    alpha(0.25)
    hold on;
    hl = line(x,y,'LineWidth',2);
    grid on
    set(hp, 'edgecolor', 'none');
    set(hl, 'color', 'k');
    xlabel('Depth (microns)','FontSize',12)
    ylabel('Autocorrelation Zero Lag','FontSize',12)
    set(gca,'Fontsize',20);



cd(od)

end



%%% SUBFUNCTIONS %%%

function [XcorrMat,lags,sig1n,sig2n] = XcorrFunction(rate1,rate2,maxlag,normMethod,spf)

% normMethod 1 is the whole matrix is normalized, normMethod 2 is each
% z-layer is normalized

Nrows = size(rate1,1);    
    
if normMethod == 1
    % normalize the entire matrix, 
    sig1n = (rate1 - mean(rate1,'all'))/std(rate1,[],'all');
    sig2n = (rate2 - mean(rate2,'all'))/std(rate2,[],'all');
end

if normMethod == 3
    % either take the avg across all layers and detrend that OR 
    % make into one continuous vector and detrend the whole thing

    % % average method
    % rate1avg = mean(rate1);
    % rate2avg = mean(rate2);
    % 
    % for zz=1:Nrows
    %     numFrames = size(rate1,2);
    %     frameVec = (1:numFrames);
    %     timeVec = (frameVec*spf)/60;
    %     sigma = 1;
    % 
    %     [sig1int(zz,:)] = interpolateAreaOrMyosin(rate1avg,sigma,timeVec,numFrames,normMethod);
    %     [sig2int(zz,:)] = interpolateAreaOrMyosin(rate2avg,sigma,timeVec,numFrames,normMethod);
    % end

    %%% continuous vector method
    rate1transpose = rate1';
    rate1vec = reshape(rate1transpose,1,[]);

    rate2transpose = rate2';
    rate2vec = reshape(rate2transpose,1,[]);

    numFrames = size(rate1vec,2);
    frameVec = (1:numFrames);
    timeVec = (frameVec*spf)/60;
    sigma = 1;

    [sig1intVec] = interpolateAreaOrMyosin(rate1vec,sigma,timeVec,numFrames,normMethod);
    [sig2intVec] = interpolateAreaOrMyosin(rate2vec,sigma,timeVec,numFrames,normMethod);

    sig1int = reshape(sig1intVec,[],size(rate1,1));
    sig1n = sig1int';

    sig2int = reshape(sig2intVec,[],size(rate2,1));
    sig2n = sig2int';

    % normalize the entire matrix, 
    %sig1n = (sig1int - mean(sig1int,'all'))/std(sig1int,[],'all');
    %sig2n = (sig2int - mean(sig2int,'all'))/std(sig2int,[],'all');


    % normalize the entire matrix, THEN detrend
    %sig1n = (rate1 - mean(rate1,'all'))/std(rate1,[],'all');
    %sig2n = (rate2 - mean(rate2,'all'))/std(rate2,[],'all');

    % for zz=1:Nrows
    %     numFrames = size(sig1norm,2);
    %     frameVec = (1:numFrames);
    %     timeVec = (frameVec*spf)/60;
    %     sigma = 1;
    % 
    %     [sig1n(zz,:)] = interpolateAreaOrMyosin(sig1norm(zz,:),sigma,timeVec,numFrames,normMethod);
    %     [sig2n(zz,:)] = interpolateAreaOrMyosin(sig2norm(zz,:),sigma,timeVec,numFrames,normMethod);
    % end

end

XcorrMat = NaN(Nrows,2*maxlag+1);
for z=1:Nrows
    
    if normMethod == 2
        sig1n(z,:) =(rate1(z,:) - mean(rate1(z,:),'all'))/std(rate1(z,:),[],'all');
        sig2n(z,:) =(rate2(z,:) - mean(rate2(z,:),'all'))/std(rate2(z,:),[],'all');
    end

    if normMethod == 4

        numFrames = size(rate1,2);
        frameVec = (1:numFrames);
        timeVec = (frameVec*spf)/60;
        sigma = 1;

        [sig1n(z,:)] = interpolateAreaOrMyosin(rate1(z,:),sigma,timeVec,numFrames,normMethod);
        [sig2n(z,:)] = interpolateAreaOrMyosin(rate2(z,:),sigma,timeVec,numFrames,normMethod);
    end
        

    [XcorrResultLayer,lags] = xcorr(sig1n(z,:),sig2n(z,:),maxlag,'normalized'); % check both 'unbiased', 'biased' and 'coeff'
    
    XcorrMat(z,:)= XcorrResultLayer;

end


end





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


function [outputArray] = interpolateAreaOrMyosin(inputArray,sigma,timeVec,numFrames,normMethod)

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

if normMethod == 3
    outputArray = ArrayVecDetrend;
end

if normMethod == 4
    % normalize by subtracting the mean and dividing by the standard deviation
    outputArray = (ArrayVecDetrend-mean(ArrayVecDetrend))/(std(ArrayVecDetrend));
end

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

