function [] = Xcorr_InterfaceMyoVsLength_LS(data,MovieNums,showSamples)
%Perform cross-correlation between the rates of junctional Myosin and
%junctional length. In this implementation for each junction we find a
%common time window in which it is present at all depths so that the
%correlation at each depth comes from the same number of time points.


% parameters 
shiftSec = 60;  % delta T for computing the rate
maxlagSec = 120; % cross-correlation maximum time lag
XYsf = 0.1;  % XY scale factor
Zsf = 0.2635; % z scale factor
fs = 10; % figure font size

zvec = 1:60; % z-layers to process
zvecMic = zvec*Zsf; % z vector in microns
Nlayers = length(zvec); % number of layers

spf = data(1).SecPerFrame;
shiftF = round(shiftSec/spf);
maxlagF = round(maxlagSec/spf);
lagvec = -maxlagF:maxlagF;
lagvecSec = lagvec*spf;
Nlags = numel(lagvec);



xcorr_resultsMat  = [];
meanAngleMat    = [];

Nmovies = length(MovieNums);
for mn=1:Nmovies
    

    cd(data(MovieNums(mn)).Source)
    cd ..

    loaddata = load('AngleShift.mat');
    AngleShift = loaddata.AngleShift;

    loadStruct = load('InterfaceBasedIntensityData.mat').InterfaceBasedIntensityData;
    InterfaceMyo = loadStruct.InterfaceIntensity_Myo(:,:,zvec);
%     Vertex1Myo   = loadStruct.Vertex1Intensity_Myo(:,:,zvec);
%     Vertex2Myo   = loadStruct.Vertex2Intensity_Myo(:,:,zvec);
%     InterfaceMyo = (Vertex1Myo + Vertex2Myo)/2;

    loadStruct = load('trackingMatrixZT_NodeFix.mat');
    Tmatrix = loadStruct.trackingMatrixZT_NodeFix(:,:,zvec);
    LengthMat = Tmatrix(:,1:8:end,:)*XYsf;
    AngleMat = Tmatrix(:,2:8:end,:)-AngleShift;
    idx = AngleMat < - 90;
    AngleMat(idx) = AngleMat(idx) + 180;
    idx = AngleMat > 90;
    AngleMat(idx) = AngleMat(idx) - 180;
    AngleMat = abs(AngleMat);

    Nints = size(LengthMat,1);

    for int = 1:Nints
        
        % get interface data and squeeze into 2D matrices where time is in
        % the column direction and depth is in the row direction.
        IntMyo = squeeze(InterfaceMyo(int,:,:))';
        IntAng = squeeze(AngleMat(int,:,:))';
        IntLen = squeeze(LengthMat(int,:,:))';

        % find a time vector in which the junction is present at all
        % depths, if it is empty of has a duration of less than 3 minutes
        % skip to next junction
        [~, tvecF] = FiniteSequenceFromNans(mean(IntMyo,1),1);
        if isempty(tvecF) || length(tvecF)*spf/60 <3
            continue;
        end

        IntMyo = IntMyo(:,tvecF);
        IntAng = IntAng(:,tvecF);
        IntLen = IntLen(:,tvecF);
        
        
        % compute rates of junction length change and junction myosin
        % intensity
        [rateMyo] = ShiftedRate( IntMyo, shiftF);
        [rateLen] = ShiftedRate( IntLen, shiftF);
        
        % perform the cross-correlation of the rates and store into a layer
        % of a results matrix.
        % normMethod 1 is the whole matrix is normalized, normMethod 2 is each
        % z-layer is normalized
        normMethod = 1;
        [xcorr_result] = XcorrFunction(rateMyo,rateLen,maxlagF,normMethod);
        xcorr_resultsMat = cat(3,xcorr_resultsMat,xcorr_result);
        
        % store time-averaged junction angles to distinguish vertical from
        % horizontal interfaces
        meanAngleVec = mean(IntAng,2,'omitnan'); 
        meanAngleMat = cat(3,meanAngleMat,meanAngleVec);

        if showSamples
            
            subplot(2,3,1)
            imagesc(tvecF*spf/60,zvecMic,IntLen);colorbar
            title('Length')
            
            subplot(2,3,2)
            imagesc(tvecF*spf/60,zvecMic,IntMyo);colorbar
            title('Myo. Junc. Int.')
            
            subplot(2,3,3)
            imagesc(lagvecSec,zvecMic,xcorr_result);colorbar
            title('Xcorr')
            
            subplot(2,3,4)
            imagesc(tvecF*spf/60,zvecMic,rateLen);colorbar
            title('Length Rate')
            
            subplot(2,3,5)
            imagesc(tvecF*spf/60,zvecMic,rateMyo);colorbar
            title('Int. Rate')
            pause
            
        end
    end % interface  
end % movies


meanXcorrL  = NaN(Nlayers,Nlags);
for z=zvec
    Az    = squeeze(meanAngleMat(z,:,:));

    XMLz = squeeze(xcorr_resultsMat(z,:,:))';
    idxA = Az < 20;
    
    meanXcorrL(z,:) = mean(XMLz(idxA,:));
end

clims = [-0.25 0.25];

subplot(2,3,1)
imagesc(lagvecSec,zvecMic,meanXcorrL)
xlabel('Lag (sec)','FontSize',fs)
ylabel('Depth (microns)','FontSize',fs)
title('Vertical: Junctional Myosin Change Vs Length Change','FontSize',fs)
colorbar
caxis(clims)
set(gca,'FontSize',fs)



subplot(2,3,4)
apical_zvec = 1:20;
lateral_zvec = 21:40;
basal_zvec  = 41:60;
plot(lagvecSec,mean(meanXcorrL(apical_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrL(lateral_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrL(basal_zvec,:),1),'LineWidth',2);hold on
xlabel('Lag (sec)','FontSize',fs)
ylabel('Cross-Correlation coeff.','FontSize',fs)
legend('Apical','Lateral','Basal')
set(gca,'FontSize',fs)
grid on
set(gca,'Ylim',clims)



meanXcorrL  = NaN(Nlayers,Nlags);

for z=zvec
    Az    = squeeze(meanAngleMat(z,:,:));

    XMLz = squeeze(xcorr_resultsMat(z,:,:))';
    idxA = Az < 60 & Az >=30;
    
    meanXcorrL(z,:) = mean(XMLz(idxA,:));
end

subplot(2,3,2)
imagesc(lagvecSec,zvecMic,meanXcorrL)
xlabel('Lag (sec)','FontSize',fs)
ylabel('Depth (microns)','FontSize',fs)
title('Diagonal: Junctional Myosin Change Vs Length Change','FontSize',fs)
colorbar
caxis(clims)
set(gca,'FontSize',fs)

subplot(2,3,5)
plot(lagvecSec,mean(meanXcorrL(apical_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrL(lateral_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrL(basal_zvec,:),1),'LineWidth',2);hold on
xlabel('Lag (sec)','FontSize',fs)
ylabel('Cross-Correlation coeff.','FontSize',fs)
legend('Apical','Lateral','Basal')
grid on
set(gca,'Ylim',clims)
set(gca,'FontSize',fs)



meanXcorrL  = NaN(Nlayers,Nlags);
meanXcorrA  = NaN(Nlayers,Nlags);
for z=zvec
    Az    = squeeze(meanAngleMat(z,:,:));
    XMLz = squeeze(xcorr_resultsMat(z,:,:))';
    idxA =  Az >=70;
    
    meanXcorrL(z,:) = mean(XMLz(idxA,:));
    

end

subplot(2,3,3)
imagesc(lagvecSec,zvecMic,meanXcorrL)
xlabel('Lag (sec)','FontSize',fs)
ylabel('Depth (microns)','FontSize',fs)
title('Horizontal: Junctional Myosin Change Vs Length Change','FontSize',fs)
colorbar
caxis(clims)
set(gca,'FontSize',fs)

subplot(2,3,6)
plot(lagvecSec,mean(meanXcorrL(apical_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrL(lateral_zvec,:),1),'LineWidth',2);hold on
plot(lagvecSec,mean(meanXcorrL(basal_zvec,:),1),'LineWidth',2);hold on
xlabel('Lag (sec)','FontSize',fs)
ylabel('Cross-Correlation coeff.','FontSize',fs)
legend('Apical','Lateral','Basal')
grid on
set(gca,'Ylim',clims)
set(gca,'FontSize',fs)



end % function


function [XcorrMat] = XcorrFunction(rate1,rate2,maxlag,normMethod)

% normMethod 1 is the whole matrix is normalized, normMethod 2 is each
% z-layer is normalized

Nrows = size(rate1,1);    
    
if normMethod == 1
    % normalize the entire matrix, 
    sig1n = (rate1 - mean(rate1,'all'))/std(rate1,[],'all');
    sig2n = (rate2 - mean(rate2,'all'))/std(rate2,[],'all');
end

XcorrMat = NaN(Nrows,2*maxlag+1);
for z=1:Nrows
    
    if normMethod == 2
        sig1n(z,:) =(rate1(z,:) - mean(rate1(z,:),'all'))/std(rate1(z,:),[],'all');
        sig2n(z,:) =(rate2(z,:) - mean(rate2(z,:),'all'))/std(rate2(z,:),[],'all');
    end
        

    [XcorrResultLayer] = xcorr(sig1n(z,:),sig2n(z,:),maxlag,'coeff'); % check both 'unbiased', 'biased' and 'coeff'
    
    XcorrMat(z,:)= XcorrResultLayer;

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

