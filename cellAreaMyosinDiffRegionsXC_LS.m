function [] = cellAreaMyosinDiffRegionsXC_LS(data,MovieNum,showSamples)

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
zvec = 1:65;
step = 1;
Alayer = 1;
Llayer = 30;%32;
Blayer = 60;%65;

Zsf = 0.2635;
zvecMic = zvec*Zsf;
XYsf = 0.1;

fs = 13;

% set lag for cross correlation
maxLag = 20;
xcorr_resultsMat  = [];

Nlags = (maxLag*2)+1;
avgXCAllMovies = [];

%shiftSec = 60;
shiftSec = 120;  % delta T for computing the rate
spf = data(1).SecPerFrame;
shiftF = round(shiftSec/spf);

xcorr_resultsMat_Int  = [];
meanAngleMat    = [];

% loop through the movies
Nmovies = numel(MovieNum);
for mn=1:Nmovies

    % get seconds per frame for the movie
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
    myoArrayI = loadMyo.CellBasedIntensityData.InteriorIntensity_Myo(:,:,zvec);
    myoArrayB = loadMyo.CellBasedIntensityData.BoundaryIntensity_Myo(:,:,zvec);
    %myoArray = loadMyo.CellBasedIntensityData.JunctionIntensity_Myo(:,:,zvec);
    %myoArray = loadMyo.CellBasedIntensityData.VertexIntensity_Myo(:,:,zvec);

    myoArray = myoArrayI + myoArrayB;



    % load myosin intensity array (interface intensity) Myosin in channel 1
    loadStruct = load('InterfaceBasedIntensityData.mat').InterfaceBasedIntensityData;
    myoArrayInt = loadStruct.InterfaceIntensity_Myo(:,:,zvec);

    % get the number of cells in the array
    Ncells = size(areaArray,1);

    % get number of layers in the array
    Nlayers = size(areaArray,3);

    count = 1;
    cmap = turbo(Nlayers);

    loaddata = load('AngleShift.mat');
    AngleShift = loaddata.AngleShift;

    loadStruct = load('trackingMatrixZT_NodeFix.mat');
    Tmatrix = loadStruct.trackingMatrixZT_NodeFix(:,:,zvec);
    LengthMat = Tmatrix(:,1:8:end,:)*XYsf;
    AngleMat = Tmatrix(:,2:8:end,:)-AngleShift;
    idx = AngleMat < - 90;
    AngleMat(idx) = AngleMat(idx) + 180;
    idx = AngleMat > 90;
    AngleMat(idx) = AngleMat(idx) - 180;
    AngleMat = abs(AngleMat);

    % loop over the cells in the movie
    for c = 1:Ncells

        % get the cell ID (inside the loop)
        CellID = cellNumbers(c,:);

        % find a time vector in which the cell is present at all
        % depths, if it is empty of has a duration of less than 3 minutes
        % skip to next cell

        cellAreaArray = squeeze(areaArray(c,:,:))'; % rows are layers, columns are time points
        cellMyoArray = squeeze(myoArray(c,:,:))';

        [~, tvecF] = FiniteSequenceFromNans(mean(cellAreaArray,1));%,1);
        if isempty(tvecF) || length(tvecF)*spf/60 <3
            continue;
        end


        % get the interface myosin (and separate into vertical / horizontal)

        % find rows in tracking matrix that contain the cell ID
        %CID1 = Tmatrix(:,7:8:end,5);
       % CID2 = Tmatrix(:,8:8:end,5);
        
        CID1 = Tmatrix(:,7:8:end,:);
        CID2 = Tmatrix(:,8:8:end,:);
        
        zeromat = CID1 == 0;
        CID1(zeromat) = NaN;
        CID2(zeromat) = NaN;
        
        CID1iscell = CID1 == CellID;
        CID2iscell = CID2 == CellID;
        cellrows = any(CID1iscell | CID2iscell, 2);
        cellrows = sum(cellrows,3);
        %cellrows = cellrows/(max(cellrows));

        cellMyoArrayInt = myoArrayInt(cellrows>0,:,:);
        AngleMatInt = AngleMat(cellrows>0,:,:);
        %cellAreaArrayInt = cellAreaArray;

         % get the cell and myosin arrays for all layers and the time vector
        cellAreaArray = cellAreaArray(:,tvecF);
        cellMyoArray = cellMyoArray(:,tvecF);
        cellMyoArrayInt = cellMyoArrayInt(:,tvecF,:);
        AngleMatInt = AngleMatInt(:,tvecF,:);

         % compute rates of cell area change and medial myosin intensity
        [cellAreaRate] = ShiftedRate( cellAreaArray, shiftF);
        [cellMyoRate] = ShiftedRate( cellMyoArray, shiftF);

        %cellAreaRate = cellAreaRate(:,71:150);
        %cellMyoRate = cellMyoRate(:,71:150);
        

        % perform the cross-correlation of the rates and store into a layer
        % of a results matrix.
        % normMethod 1 is the whole matrix is normalized, normMethod 2 is each
        % z-layer is normalized
        normMethod = 1;
        %[xcorr_result,lags,areaNorm,myoNorm] = XcorrFunction(cellAreaArray,cellMyoArray,maxLag,normMethod,spf);
        [xcorr_result,lags,areaNorm,myoNorm] = XcorrFunction(cellAreaRate,cellMyoRate,maxLag,normMethod,spf);
        xcorr_resultsMat = cat(3,xcorr_resultsMat,xcorr_result);

        if MovieNum(mn) == 3 && c == 79
            xcorr_CellEx = xcorr_result;
        end

        lagvec = -maxLag:maxLag;
        lagvecSec = lagvec*spf;

        % loop through interfaces and store the cross correlation results

        Nints = size(cellMyoArrayInt,1);
        for i = 1:Nints

            intMyo = squeeze(cellMyoArrayInt(i,:,:))';
            intAng = squeeze(AngleMatInt(i,:,:))';

            [~, tvecF] = FiniteSequenceFromNans(mean(intMyo,1));%,1);
            if isempty(tvecF) || length(tvecF)*spf/60 <3
                continue;
            end

            intMyo = intMyo(:,tvecF);
            intAng = intAng(:,tvecF);
            intCellArea = cellAreaArray(:,tvecF);


            [intMyoRate] = ShiftedRate( intMyo, shiftF);
            [intCellAreaRate] = ShiftedRate( intCellArea, shiftF);

            % intMyoRate = intMyo;
            % intCellAreaRate = intCellArea;

            % store junction angles
            meanAngleVec = mean(intAng,2,'omitnan'); 
            meanAngleMat = cat(3,meanAngleMat,meanAngleVec);


            [xcorr_result_Int,lags,areaIntNorm,myoIntNorm] = XcorrFunction(intCellAreaRate,intMyoRate,maxLag,normMethod,spf);
            xcorr_resultsMat_Int = cat(3,xcorr_resultsMat_Int,xcorr_result_Int);

            % if showSamples
            % 
            %     subplot(3,1,1)
            %     imagesc(tvecF*spf/60,zvecMic,intCellAreaRate);colorbar
            %     title('Area')
            % 
            %     subplot(3,1,2)
            %     imagesc(tvecF*spf/60,zvecMic,intMyoRate);colorbar
            %     title('Myosin')
            % 
            %     subplot(3,1,3)
            %     imagesc(lagvecSec,zvecMic,xcorr_result_Int);colorbar
            %     title('Xcorr')
            % 
            %     %pause
            % 
            % end

        end


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


        avgXC  = NaN(Nlayers,Nlags);
        % for l = 1:step:Nlayers
        % 
        %     XMLz = squeeze(xcorr_resultsMat(l,:,:))';
        % 
        %     avgXC(l,:) = mean(XMLz,'omitnan');
        % 
        %     % plot the average cross correlation
        %     plot(lags*spf/60,avgXC(l,:),'LineWidth',2,'Color',cmap(l,:));
        %     set(gca,'FontSize',16)
        %     txt = ['Movie ',num2str(MovieNum),' (blue to red in layers)'];%' Z depth ',num2str(zSpacing*(j-1)),' microns'];
        %     title(txt)
        %     xlabel('Lag (minutes)')
        %     ylabel('Correlation')
        %     grid on
        %     %ylim([-40 150])
        %     %pause
        %     hold on
        %     set(gca,'Fontsize',20)
        %     %legend('0 microns','4 microns','8 microns','12 microns','16 microns')
        % 
        %     count = count + step;
        % 
        % end


    end

    %avgXCAllMovies = cat(3,avgXCAllMovies,avgXC);

    zeroLagMovie = squeeze(xcorr_resultsMat(:,maxLag+1,:));
    zeroLagXC{:,mn} = zeroLagMovie;

    % reinitialize arrays for averaging
    avgXC = [];
    avgLags = [];
    count = 1;
    %xcorr_resultsMat = [];

end

%avgXCAllMovies = mean(avgXCAllMovies,3,'omitnan');

meanXcorrI  = NaN(Nlayers,Nlags);
for z=zvec
    XMLz = squeeze(xcorr_resultsMat(z,:,:))';
    
    meanXcorrI(z,:) = mean(XMLz,'omitnan');
end


lagvecSec = lags*spf;

figure
subplot(1,2,1)
imagesc(lagvecSec,zvecMic,meanXcorrI)
xlabel('Lag (sec)')
ylabel('Depth (microns)')
%title('Cross Sectional Area vs Medial Myosin')
title('Area vs Interior Myosin Rates')
colorbar
caxis([-0.4 0.1])
set(gca,'FontSize',20)

subplot(1,2,2)
%top_zvec = 1;
%apical_zvec = 1;%round(5/Zsf);
%lateral_zvec = round(8/Zsf);%round(10/Zsf);
%basal_zvec  = 60;%round(16/Zsf);%round(15/Zsf);


 apical_zvec = 1:20;
 lateral_zvec = 21:40;
 basal_zvec  = 41:60;
%plot(lagvecSec,mean(meanXcorrI(top_zvec,:),1),'LineWidth',2);hold on
%plot(lagvecSec,mean(meanXcorrI(apical_zvec,:),1),'LineWidth',2);hold on
%plot(lagvecSec,mean(meanXcorrI(lateral_zvec,:),1),'LineWidth',2);hold on
%plot(lagvecSec,mean(meanXcorrI(basal_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI,1),'LineWidth',2)
xlabel('Lag (sec)')
ylabel('Cross-Correlation coeff.')
%legend('0 microns','5 microns','10 microns','15 microns')
%legend('0 microns','8 microns','16 microns')
legend('Apical','Lateral','Basal')
set(gca,'FontSize',20)
grid on
set(gca,'Ylim',[-0.4 0.1])

hold on
plot(lagvecSec,mean(xcorr_CellEx,1),'LineWidth',2)
legend('All Cells','Example Cell')


% interface myosin (vertical)
meanXcorrI  = NaN(Nlayers,Nlags);
for z=zvec
    Az    = squeeze(meanAngleMat(z,:,:));

    XMLz = squeeze(xcorr_resultsMat_Int(z,:,:))';
    idxA = Az < 20;
    
    meanXcorrI(z,:) = mean(XMLz(idxA,:),'omitnan');
end

clims = [-0.4 0.1];

figure
subplot(2,3,1)
imagesc(lagvecSec,zvecMic,meanXcorrI)
xlabel('Lag (sec)')
ylabel('Depth (microns)')
title('Vertical: Junctional Myosin Vs Area Rates')
colorbar
caxis(clims)
set(gca,'FontSize',fs)



subplot(2,3,4)
top_zvec = 1;
apical_zvec = round(5/Zsf);
lateral_zvec = round(10/Zsf);
basal_zvec  = round(15/Zsf);
plot(lagvecSec,mean(meanXcorrI(top_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(apical_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(lateral_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(basal_zvec,:),1),'LineWidth',2);hold on
xlabel('Lag (sec)')
ylabel('Cross-Correlation coeff.')
legend('0 microns','5 microns','10 microns','15 microns')
grid on
set(gca,'Ylim',clims)
set(gca,'FontSize',fs)


% interface myosin (transverse)
meanXcorrI  = NaN(Nlayers,Nlags);
for z=zvec
    Az    = squeeze(meanAngleMat(z,:,:));

    XMLz = squeeze(xcorr_resultsMat_Int(z,:,:))';
    idxA = Az < 60 & Az >= 30;
    
    meanXcorrI(z,:) = mean(XMLz(idxA,:),'omitnan');
end



subplot(2,3,2)
imagesc(lagvecSec,zvecMic,meanXcorrI)
xlabel('Lag (sec)')
ylabel('Depth (microns)')
title('Transverse: Junctional Myosin Vs Area Rates')
colorbar
caxis(clims)
set(gca,'FontSize',fs)



subplot(2,3,5)
top_zvec = 1;
apical_zvec = round(5/Zsf);
lateral_zvec = round(10/Zsf);
basal_zvec  = round(15/Zsf);
plot(lagvecSec,mean(meanXcorrI(top_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(apical_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(lateral_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(basal_zvec,:),1),'LineWidth',2);hold on
xlabel('Lag (sec)')
ylabel('Cross-Correlation coeff.')
legend('0 microns','5 microns','10 microns','15 microns')
grid on
set(gca,'Ylim',clims)
set(gca,'FontSize',fs)

% interface myosin (horizontal)
meanXcorrI  = NaN(Nlayers,Nlags);
for z=zvec
    Az    = squeeze(meanAngleMat(z,:,:));

    XMLz = squeeze(xcorr_resultsMat_Int(z,:,:))';
    idxA = Az >= 70;
    
    meanXcorrI(z,:) = mean(XMLz(idxA,:),'omitnan');
end



subplot(2,3,3)
imagesc(lagvecSec,zvecMic,meanXcorrI)
xlabel('Lag (sec)')
ylabel('Depth (microns)')
title('Horizontal: Junctional Myosin Vs Area Rates')
colorbar
caxis(clims)
set(gca,'FontSize',fs)



subplot(2,3,6)
top_zvec = 1;
apical_zvec = round(5/Zsf);
lateral_zvec = round(10/Zsf);
basal_zvec  = round(15/Zsf);
plot(lagvecSec,mean(meanXcorrI(top_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(apical_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(lateral_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrI(basal_zvec,:),1),'LineWidth',2);hold on
xlabel('Lag (sec)')
ylabel('Cross-Correlation coeff.')
legend('0 microns','5 microns','10 microns','15 microns')
grid on
set(gca,'Ylim',clims)
set(gca,'FontSize',fs)


end



%%% SUBFUNCTIONS %%%

function [XcorrMat,lags,sig1n,sig2n] = XcorrFunction(rate1,rate2,maxlag,normMethod,spf)

% normMethod 1 is the whole matrix is normalized, normMethod 2 is each
% z-layer is normalized

Nrows = size(rate1,1);    
    
if normMethod == 1
    % normalize the entire matrix, 
    sig1n = (rate1 - mean(rate1,'all','omitnan'))/std(rate1,[],'all','omitnan');
    sig2n = (rate2 - mean(rate2,'all','omitnan'))/std(rate2,[],'all','omitnan');
end

XcorrMat = NaN(Nrows,2*maxlag+1);
for z=1:Nrows
    
    if normMethod == 2
        sig1n(z,:) =(rate1(z,:) - mean(rate1(z,:),'all'))/std(rate1(z,:),[],'all');
        sig2n(z,:) =(rate2(z,:) - mean(rate2(z,:),'all'))/std(rate2(z,:),[],'all');
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






