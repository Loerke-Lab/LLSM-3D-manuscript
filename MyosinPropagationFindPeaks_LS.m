function [ExpandMatAll,ContractMatAll] = MyosinPropagationFindPeaks_LS(data,MovieNums,showEgs)
%This function analyzes myosin oscillation z-propagation direction and speed.
%
% INPUTS:
%           data: structure that contains Source, ImageFileList, and
%                   SecPerFrame fields for each movie
%           MovieNums: vector of movie indices (of data) to include in
%                      analysis
%           showEgs: to see visualization set to 1, else set to 0.
%
%
% OUTPUTS:
%           ExpandMatAll:  MxNx2 Matrix that contains Expansion peak data.
%           M is the number of z-layers, N is the number of peaks found and
%           the 2 layers correspond to Peak Myosin Rate (1st layer), and Peak
%           Myosin time (2nd layer).


close all;

zvec = 1:60;
Nz = numel(zvec);

% length and time scale factors
XYsf = 0.1;
Zsf  = 0.2635;
spf = data(MovieNums(1)).SecPerFrame;

% rate parameter (delta T)
DeltaTframes = round(30/spf);
DeltaTminutes = DeltaTframes*spf/60;

% minimum duration in minutes for the cell
minDurationMin = 4;
minDurationFrames =round(minDurationMin*60/spf);

% initialize results matrices
ExpandMatAll = [];
ContractMatAll = [];

CellNumber = 1;
Nmovies = numel(MovieNums);
for mn=1:Nmovies

    % move to data directory
    cd(data(MovieNums(mn)).Source)
    cd ..

    % load myosin intensity array (interior - M for medial) Myosin in channel 1
    loadMyo = load('CellBasedIntensityData.mat');
    myoArrayI = loadMyo.CellBasedIntensityData.InteriorIntensity_Myo(:,:,zvec);
    myoArrayB = loadMyo.CellBasedIntensityData.BoundaryIntensity_Myo(:,:,zvec);
    %myoArray = loadMyo.CellBasedIntensityData.JunctionIntensity_Myo(:,:,zvec);

    myoArray = myoArrayI + myoArrayB;


    % % load cell area data (TO TEST CODE)
    % loadGdata = load('GeometryData.mat');
    % Area3D = loadGdata.Area(:,:,zvec)*XYsf^2;
    % myoArray = Area3D;


    Ncells = size(myoArray,1); % number of cells to loop over

    for icell = 1:Ncells

        MyoZvT = squeeze(myoArray(icell,:,:))'; % squeeze cell area so we have depth (rows) versus time (columns)

        % find the longest finite sequence in which all z-layers are
        % present
        [~,frameVec] = FiniteSequenceFromNans(mean(MyoZvT));%,1);
        
        % if time course is too short skip over it
        if numel(frameVec)< minDurationFrames
            continue;
        end

        % crop cell area matrix to frame vector and convert frame vector to
        % time vector in minutes
        MyoZvT = MyoZvT(:,frameVec);
        tvecMin = frameVec*spf/60;


        % median filtering is used to smooth over transient detection errors
        MyoZvTmed = MyoZvT;%medfilt2(AreaZvT,[5 5],'symmetric');
        
        % take the rate of change
        MyoZvT2 = MyoZvTmed(:,(DeltaTframes+1):end);
        MyoZvT1 = MyoZvTmed(:,1:(end-DeltaTframes));
        MyoRateZvT = (MyoZvT2 - MyoZvT1)/DeltaTminutes;
        
        % filter the rate so that findpeaks finds true peaks and not noise
        [MyoRateZvT] = filterImage3DpaddedEdges(MyoRateZvT, 'Gauss', [3,2,0]);

        % % normalize myosin
        % normMethod = 1;
        % if normMethod == 1
        %     % normalize the entire matrix, 
        %     MyoRateZvT = (MyoRateZvT - mean(MyoRateZvT,'all'))/std(MyoRateZvT,[],'all');
        % 
        % end

        % initialize logical Z vs. T matrices where we will 'true' the peak
        % positions for each depth.
        expandPeaksLogicalZvT = false(size(MyoRateZvT));
        cntractPeaksLogicalZvT = false(size(MyoRateZvT));
        for z=1:Nz
            [~,locs] = findpeaks(MyoRateZvT(z,:));             
            expandPeaksLogicalZvT(z,locs) = true;

            [~,locs] = findpeaks(-MyoRateZvT(z,:));           
            cntractPeaksLogicalZvT(z,locs) = true;
        end

        % now label each 'crest', first dilate them to connect any that are
        % very close, we will undo the dilation afterwards
        exp_dil=imdilate(expandPeaksLogicalZvT,strel('disk',1));
        cnt_dil=imdilate(cntractPeaksLogicalZvT,strel('disk',1));
        Lexp = bwlabel(exp_dil);
        Lcnt = bwlabel(cnt_dil);
        % now zero all the pixels that are not peaks (i.e. undo dilation pixels)
        Lexp(~expandPeaksLogicalZvT) = 0;
        Lcnt(~cntractPeaksLogicalZvT) = 0;

        % get number of labels for both contractions and expansions and for
        % each loop over each one storing the rates (not propagation rate but myosin oscillation rate) and time points along
        % the labeled peak
        Nexp = max(Lexp,[],'all');
        ExpandMatIndividual = NaN(Nz,Nexp,2);
        for ii=1:Nexp
            bw = Lexp == ii;
            ind = find(bw);
            [row,col] = ind2sub(size(bw),ind);
            ratevec = MyoRateZvT(ind);
            [B,I] = sort(row);
            if any(diff(B)==0) % branch exists, so skip to next
                continue;
            end
            RateVec = NaN(Nz,1);
            RateVec(B) = ratevec(I);
            TimeVec = NaN(Nz,1);
            TimeVec(B) = col(I);
            ExpandMatIndividual(:,ii,1) = RateVec;
            ExpandMatIndividual(:,ii,2) = TimeVec;
        end


        Ncnt = max(Lcnt,[],'all');
        ContractMatIndividual = NaN(Nz,Ncnt,2);
        for ii=1:Ncnt
            bw = Lcnt == ii;
            ind = find(bw);
            [row,col] = ind2sub(size(bw),ind);
            ratevec = MyoRateZvT(ind);
            [B,I] = sort(row);
            if any(diff(B)==0) % branch exists, so skip to next
                continue;
            end
            RateVec = NaN(Nz,1);
            RateVec(B) = ratevec(I);
            TimeVec = NaN(Nz,1);
            TimeVec(B) = double(col(I));
            ContractMatIndividual(:,ii,1) = RateVec;
            ContractMatIndividual(:,ii,2) = TimeVec;
        end

        ExpandMatAll = [ExpandMatAll,ExpandMatIndividual];
        ContractMatAll = [ContractMatAll,ContractMatIndividual];
        CellNumber = CellNumber + 1;


        if showEgs

            subplot(1,3,1)
            imagesc(tvecMin,zvec*Zsf,MyoRateZvT);colorbar

            subplot(1,3,2)
            ADFe = MyoRateZvT;
            ADFe(Lexp>0) = NaN;
            imagesc(tvecMin,zvec*Zsf,ADFe)
            colorbar;
            title('Expansion Peaks')
            
            
            subplot(1,3,3)
            imagesc(tvecMin,zvec*Zsf,Lexp)
            title('Expansion Peaks Labeled')
             
            pause;
        end


    end % cell
    
    
end % movie
%close all

figure (1)
%edges = -120:20:120;
%centers = -110:20:110;

% edges = -110:20:110;
% centers = -100:20:100;

% edges = -100:10:100;
% centers = -95:10:95;


%%% Added 4/26/23 KB %%%
% centers = -60:10:60;
% centers = -120:20:120;
% centers = -(spf*9):spf:(spf*9);

%centers = -90:15:90;
%centers = -100:10:100;

% d = diff(centers)/2;
% edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];

sizeHist = 22;
sizeJump = 38; % 10 microns

% sizeHist = 41;
% sizeJump = 19; % 5 microns

% sizeHist = 30;
% sizeJump = 30; % 8 microns

% sizeHist = 37;
% sizeJump = 23; % 6 microns

% sizeHist = 45;
% sizeJump = 15; % 4 microns

% %%% Added 5/1/23 KB %%%
% Nedges = length(edges)-1;
% ApicalExpHist = NaN(sizeHist,Nedges);
% 
% for ii=1:sizeHist
%     DiffBasalExp = spf*(ExpandMatAll(ii+sizeJump,:,2)-ExpandMatAll(ii,:,2)); % expansion (myosin increase)
%     %DiffBasalExp = spf*(ContractMatAll(ii+sizeJump,:,2)-ContractMatAll(ii,:,2)); % contraction (myosin decrease)
%     DiffBasalExp(isnan(DiffBasalExp)) = [];
%     %DiffBasalExp(DiffBasalExp==0) = [];
%     [N,~]= histcounts(DiffBasalExp,edges);%,'Normalization', 'probability');
% 
%     ApicalExpHist(ii,:) = N;
% 
%     % if ii==1 | ii==11 | ii==22
%     %     figure(3)
%     %     histogram(DiffBasalExp,edges)
%     %     hold on
%     % end
% end
%%%

centers = -(spf*12):spf:(spf*12);
d = diff(centers)/2;
edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];


Nedges = length(edges)-1;
ApicalExpHist = NaN(22,Nedges);

posValuesA = [];
posValuesB = [];
negValuesA = [];
negValuesB = [];

for ii=1:22
    DiffBasalExp = spf*(ExpandMatAll(ii+38,:,2)-ExpandMatAll(ii,:,2));
    DiffBasalExp(isnan(DiffBasalExp)) = [];
    %DiffBasalExp(DiffBasalExp==0) = []; % added for additional analysis, can remove
    [N,~,bin]= histcounts(DiffBasalExp,edges,'Normalization', 'probability');

    ApicalExpHist(ii,:) = N;

    % save positive and negative values for innerquartile range later (Added KL 5/8/2024)
    if ii >= 1 && ii <=5
        positiveValues = DiffBasalExp(DiffBasalExp>0);
        posValuesA = horzcat(posValuesA,positiveValues);

        negativeValues = DiffBasalExp(DiffBasalExp<0);
        negValuesA = horzcat(negValuesA,negativeValues);

        positiveValues = [];
        negativeValues = [];
    elseif ii >= 18 && ii <= 22
        positiveValues = DiffBasalExp(DiffBasalExp>0);
        posValuesB = horzcat(posValuesB,positiveValues);

        negativeValues = DiffBasalExp(DiffBasalExp<0);
        negValuesB = horzcat(negValuesB,negativeValues);

        positiveValues = [];
        negativeValues = [];
    end


end


apicalLayers = 1:5; % first five
basalLayers = 18:22; % last five

lateralLayers = 11;


% percent positive and negative
numColumns = numel(centers);
% get the percent negative and percent positive values
binsNegative = ApicalExpHist(:,1:round(numColumns/2)-1);
binsPositive = ApicalExpHist(:,round(numColumns/2)+1:end);
binsAll = horzcat(ApicalExpHist(:,1:round(numColumns/2)-1),ApicalExpHist(:,round(numColumns/2)+1:end));
prcPositive = ((sum(binsPositive,'all'))/(sum(binsAll,'all')))*100;
prcNegative = ((sum(binsNegative,'all'))/(sum(binsAll,'all')))*100;
prcPositive = ((sum(binsPositive,2))./(sum(binsAll,2)))*100;
prcNegative = ((sum(binsNegative,2))./(sum(binsAll,2)))*100;
% apical section
prcPositiveA = ((sum(binsPositive(apicalLayers,:),'all'))/(sum(binsAll(apicalLayers,:),'all')))*100;
prcNegativeA = ((sum(binsNegative(apicalLayers,:),'all'))/(sum(binsAll(apicalLayers,:),'all')))*100;
% basal section
prcPositiveB = ((sum(binsPositive(basalLayers,:),'all'))/(sum(binsAll(basalLayers,:),'all')))*100;
prcNegativeB = ((sum(binsNegative(basalLayers,:),'all'))/(sum(binsAll(basalLayers,:),'all')))*100;

% lateral section
prcPositiveL = ((sum(binsPositive(lateralLayers,:),'all'))/(sum(binsAll(lateralLayers,:),'all')))*100;
prcNegativeL = ((sum(binsNegative(lateralLayers,:),'all'))/(sum(binsAll(lateralLayers,:),'all')))*100;


% find center of gravity for negative and positive
timeVecNeg = centers(:,1:round(numColumns/2)-1);
timeVecPos = centers(:,round(numColumns/2)+1:end);

%timeVecNeg = repmat(timeVecNeg,5,1);
%timeVecPos = repmat(timeVecPos,5,1);

% APICAL
binsNegativeA = mean(binsNegative(apicalLayers,:),1);
binsPositiveA = mean(binsPositive(apicalLayers,:),1);

sumPositiveA = sum(binsPositiveA,'all');
sumNegativeA = sum(binsNegativeA,'all');

sumAmpPositiveA = sum(binsPositiveA.*timeVecPos,'all');
sumAmpNegativeA = sum(binsNegativeA.*timeVecNeg,'all');

cogPositiveA = sumAmpPositiveA/sumPositiveA;
cogNegativeA = sumAmpNegativeA/sumNegativeA;

% get y value for center of gravity
idxPosA1 = find(timeVecPos < cogPositiveA,1,'last');
idxPosA2 = find(timeVecPos > cogPositiveA,1,'first');
idxNegA1 = find(timeVecNeg < cogNegativeA,1,'last');
idxNegA2 = find(timeVecNeg > cogNegativeA,1,'first');

probPosA = interp1(timeVecPos(idxPosA1:idxPosA2),binsPositiveA(idxPosA1:idxPosA2),cogPositiveA);
probNegA = interp1(timeVecNeg(idxNegA1:idxNegA2),binsNegativeA(idxNegA1:idxNegA2),cogNegativeA);

% BASAL
binsNegativeB = mean(binsNegative(basalLayers,:),1);
binsPositiveB = mean(binsPositive(basalLayers,:),1);

sumPositiveB = sum(binsPositiveB,'all');
sumNegativeB = sum(binsNegativeB,'all');

sumAmpPositiveB = sum(binsPositiveB.*timeVecPos,'all');
sumAmpNegativeB = sum(binsNegativeB.*timeVecNeg,'all');

cogPositiveB = sumAmpPositiveB/sumPositiveB;
cogNegativeB = sumAmpNegativeB/sumNegativeB;

% get y value for center of gravity
idxPosB1 = find(timeVecPos < cogPositiveB,1,'last');
idxPosB2 = find(timeVecPos > cogPositiveB,1,'first');
idxNegB1 = find(timeVecNeg < cogNegativeB,1,'last');
idxNegB2 = find(timeVecNeg > cogNegativeB,1,'first');

probPosB = interp1(timeVecPos(idxPosB1:idxPosB2),binsPositiveB(idxPosB1:idxPosB2),cogPositiveB);
probNegB = interp1(timeVecNeg(idxNegB1:idxNegB2),binsNegativeB(idxNegB1:idxNegB2),cogNegativeB);

% LATERAL
binsNegativeL = mean(binsNegative(lateralLayers,:),1);
binsPositiveL = mean(binsPositive(lateralLayers,:),1);

sumPositiveL = sum(binsPositiveL,'all');
sumNegativeL = sum(binsNegativeL,'all');

sumAmpPositiveL = sum(binsPositiveL.*timeVecPos,'all');
sumAmpNegativeL = sum(binsNegativeL.*timeVecNeg,'all');

cogPositiveL = sumAmpPositiveL/sumPositiveL;
cogNegativeL = sumAmpNegativeL/sumNegativeL;

% get y value for center of gravity
idxPosL1 = find(timeVecPos < cogPositiveL,1,'last');
idxPosL2 = find(timeVecPos > cogPositiveL,1,'first');
idxNegL1 = find(timeVecNeg < cogNegativeL,1,'last');
idxNegL2 = find(timeVecNeg > cogNegativeL,1,'first');

probPosL = interp1(timeVecPos(idxPosL1:idxPosL2),binsPositiveL(idxPosL1:idxPosL2),cogPositiveL);
probNegL = interp1(timeVecNeg(idxNegL1:idxNegL2),binsNegativeL(idxNegL1:idxNegL2),cogNegativeL);


% get the range of speeds for 50% of the data
[iqrPosA,xValsPosA] = iqr(posValuesA);
binPosA = find(fix(centers)==fix(xValsPosA(1)) | fix(centers)==fix(xValsPosA(2)));
yValsPosA = mean(ApicalExpHist(apicalLayers,binPosA));

[iqrPosB,xValsPosB] = iqr(posValuesB);
binPosB = find(fix(centers)==fix(xValsPosB(1)) | fix(centers)==fix(xValsPosB(2)));
yValsPosB = mean(ApicalExpHist(basalLayers,binPosB));

[iqrNegA,xValsNegA] = iqr(negValuesA);
binNegA = find(fix(centers)==fix(xValsNegA(1)) | fix(centers)==fix(xValsNegA(2)));
yValsNegA = mean(ApicalExpHist(apicalLayers,binNegA));

[iqrNegB,xValsNegB] = iqr(negValuesB);
binNegB = find(fix(centers)==fix(xValsNegB(1)) | fix(centers)==fix(xValsNegB(2)));
yValsNegB = mean(ApicalExpHist(basalLayers,binNegB));


% PLOTS
subplot(1,3,1)
imagesc(centers,(1:60)*0.2635,ApicalExpHist);title('Contraction from apical end')
title('\Delta z = 10 \mum')
xlabel('\Delta t (s)')
ylabel('Depth (\mum)')
colorbar
clim([0 0.14])
set(gca,'Fontsize',18)
xlim([-(spf*12) spf*12])
xticks([-80 -60 -40 -20 0 20 40 60 80])



subplot(1,3,2)
%plot(centers,mean(ApicalExpHist(apicalLayers,:),1),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
a = area(centers(1:13),mean(ApicalExpHist(apicalLayers,1:13),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%EdgeColor','none')
hold on
b = area(centers(13:25),mean(ApicalExpHist(apicalLayers,13:25),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%'EdgeColor','none')
a.FaceAlpha = 0.5;
b.FaceAlpha = 0.5;
xlabel('\Delta t (s)')
ylabel('Probability')
set(gca,'Fontsize',18)
ylim([0 0.14])
xlim([-(spf*12) spf*12])
xticks([-80 -60 -40 -20 0 20 40 60 80])
grid on
%alpha(0.5)
c = area(centers(binNegA(1):binNegA(2)),mean(ApicalExpHist(apicalLayers,binNegA(1):binNegA(2)),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%EdgeColor','none'))
d = area(centers(binPosA(1):binPosA(2)),mean(ApicalExpHist(apicalLayers,binPosA(1):binPosA(2)),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%EdgeColor','none'))
%c.FaceAlpha = 0.75;
%d.FaceAlpha = 0.75;
scatter(cogNegativeA,probNegA,75,'k','filled')
scatter(cogPositiveA,probPosA,75,'k','filled')

%plot([xValsNegA(1),xValsNegA(1)], [0,yValsNegA(1)],'k','LineWidth',1)
%plot([xValsNegA(2),xValsNegA(2)], [0,yValsNegA(2)],'k','LineWidth',1)
%plot([xValsPosA(1),xValsPosA(1)], [0,yValsPosA(1)],'k','LineWidth',1)
%plot([xValsPosA(2),xValsPosA(2)], [0,yValsPosA(2)],'k','LineWidth',1)

% subplot(1,4,3)
% %plot(centers,mean(ApicalContHist(apicalLayers,:),1),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% area(centers(1:13),mean(ApicalExpHist(lateralLayers,1:13),1),'FaceColor',[0 0.52 0.54],'LineWidth',2)%EdgeColor','none')
% hold on
% area(centers(13:25),mean(ApicalExpHist(lateralLayers,13:25),1),'FaceColor',[0.89 0.43 0],'LineWidth',2)%'EdgeColor','none')
% xlabel('\Delta t (s)')
% ylabel('Probability')
% set(gca,'Fontsize',18)
% ylim([0 0.14])
% xlim([-(spf*12) spf*12])
% xticks([-80 -60 -40 -20 0 20 40 60 80])
% grid on
% alpha(0.5)
% scatter(cogNegativeL,probNegL,75,'k','filled')
% scatter(cogPositiveL,probPosL,75,'k','filled')

subplot(1,3,3)
%plot(centers,mean(ApicalExpHist(basalLayers,:),1),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
a = area(centers(1:13),mean(ApicalExpHist(basalLayers,1:13),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%'EdgeColor','none')
hold on
b = area(centers(13:25),mean(ApicalExpHist(basalLayers,13:25),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%'EdgeColor','none')
a.FaceAlpha = 0.5;
b.FaceAlpha = 0.5;
xlabel('\Delta t (s)')
ylabel('Probability')
set(gca,'Fontsize',18)
ylim([0 0.14])
xlim([-(spf*12) spf*12])
xticks([-80 -60 -40 -20 0 20 40 60 80])
grid on
%alpha(0.5)
c = area(centers(binNegB(1):binNegB(2)),mean(ApicalExpHist(basalLayers,binNegB(1):binNegB(2)),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%EdgeColor','none'))
d = area(centers(binPosB(1):binPosB(2)),mean(ApicalExpHist(basalLayers,binPosB(1):binPosB(2)),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%EdgeColor','none'))
%c.FaceAlpha = 0.75;
%d.FaceAlpha = 0.75;
scatter(cogNegativeB,probNegB,75,'k','filled')
scatter(cogPositiveB,probPosB,75,'k','filled')

%plot([xValsNegB(1),xValsNegB(1)], [0,yValsNegB(1)],'k','LineWidth',1)
%plot([xValsNegB(2),xValsNegB(2)], [0,yValsNegB(2)],'k','LineWidth',1)
%plot([xValsPosB(1),xValsPosB(1)], [0,yValsPosB(1)],'k','LineWidth',1)
%plot([xValsPosB(2),xValsPosB(2)], [0,yValsPosB(2)],'k','LineWidth',1)


end


% subplot(1,4,1)
% %imagesc(edges(1:end-1),(1:60)*0.2635,ApicalExpHist);title('Expansion from apical end')
% %imagesc(edges,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
% imagesc(centers,(1:60)*0.2635,ApicalExpHist);title('Contraction from apical end')
% title('\Delta z = 10 \mum','FontSize',16)
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Depth (\mum)','FontSize',16)
% colorbar
% %caxis([0 0.35])
% %xlim([-100 100])
% %xlim([-49 49])
% %xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
% %xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])
% %caxis([0 0.20]);
% %yticks([0 5 10 15])
% %set(gca,'FontSize',20)
% 
% 
% %subplot(1,4,2)
% % plot(edges(1:end-1),ApicalExpHist([1 11 22],:))
% % legend('Apical 10','Lateral 10','Basal 10')
% % xlabel('\Delta t (s)','FontSize',16)
% % ylabel('Probability','FontSize',16)
% % title('+ = Apical leading')
% 
% 
% subplot(1,4,2)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalExpHist(1,:),[0.1 -50 22 0.2 50 22]);
% %estimate = fitcurveMultiGauss1D_simple(centers,ApicalExpHist(1,:),[0.1 -50 22 0.2 50 22]);
% plot(centers,ApicalExpHist(1,:),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% xlim([-100 100])
% grid on
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% subplot(1,4,3)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalExpHist(11,:),[0.1 -50 22 0.2 50 22]);
% %estimate = fitcurveMultiGauss1D_simple(centers,ApicalExpHist(11,:),[0.1 -50 22 0.2 50 22]);
% plot(centers,ApicalExpHist(11,:),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% xlim([-100 100])
% grid on
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% subplot(1,4,4)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalExpHist(22,:),[0.1 -50 22 0.2 50 22]);
% %estimate = fitcurveMultiGauss1D_simple(centers,ApicalExpHist(22,:),[0.1 -50 22 0.2 50 22]);
% plot(centers,ApicalExpHist(22,:),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% xlim([-100 100])
% grid on
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% % figure
% % plot(edges(1:end-1),ApicalExpHist(1,:),'Linewidth',2)
% % hold on
% % plot(edges(1:end-1),ApicalExpHist(11,:),'Linewidth',2)
% % plot(edges(1:end-1),ApicalExpHist(22,:),'Linewidth',2)
% % legend('0-10 microns','3-13 microns','6-16 microns')
% % title('Myosin Propagation')
% % set(gca,'Fontsize',18)
% % xlim([-100 100])
% % grid on
% % 
% % figure
% % histogram(ApicalExpHist(1,:),20)
% % hold on
% % histogram(ApicalExpHist(11,:),20)
% % histogram(ApicalExpHist(22,:),20)
% % legend('0-10 microns','3-13 microns','6-16 microns')
% % title('Myosin Propagation')
% % set(gca,'Fontsize',18)
% % xlim([-100 100])
% % grid on

%end















