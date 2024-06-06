function [ExpandMatAll,ContractMatAll] = AreaPropagationFindPeaks_LS_revised(data,MovieNums,showEgs)
%This function analyzes area oscillation z-propagation direction and speed.
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
%           the 2 layers correspond to Peak Area Rate (1st layer), and Peak
%           Area time (2nd layer).

%close all;
%figure

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
    
    % load cell area data
    loadGdata = load('GeometryData.mat');
    Area3D = loadGdata.Area(:,:,zvec)*XYsf^2;
    Ncells = size(Area3D,1); % number of cells to loop over
    for icell = 1:Ncells
        
        AreaZvT = squeeze(Area3D(icell,:,:))'; % squeeze cell area so we have depth (rows) versus time (columns)
        
        % find the longest finite sequence in which all z-layers are
        % present
        [~,frameVec] = FiniteSequenceFromNans(mean(AreaZvT));%,1);
        
        % if time course is too short skip over it
        if numel(frameVec)< minDurationFrames
            continue;
        end
        
        % crop cell area matrix to frame vector and convert frame vector to
        % time vector in minutes
        AreaZvT = AreaZvT(:,frameVec);
        tvecMin = frameVec*spf/60;
        
        % median filtering is used to smooth over transient detection errors
        AreaZvTmed = AreaZvT;%medfilt2(AreaZvT,[5 5],'symmetric');
        
        % take the rate of change
        AreaZvT2 = AreaZvTmed(:,(DeltaTframes+1):end);
        AreaZvT1 = AreaZvTmed(:,1:(end-DeltaTframes));
        AreaRateZvT = (AreaZvT2 - AreaZvT1)/DeltaTminutes;
        
        % filter the rate so that findpeaks finds true peaks and not noise
        [AreaRateZvT] = filterImage3DpaddedEdges(AreaRateZvT, 'Gauss', [3,2,0]);


        % initialize logical Z vs. T matrices where we will 'true' the peak
        % positions for each depth.
        expandPeaksLogicalZvT = false(size(AreaRateZvT));
        cntractPeaksLogicalZvT = false(size(AreaRateZvT));
        for z=1:Nz
            [~,locs] = findpeaks(AreaRateZvT(z,:));             
            expandPeaksLogicalZvT(z,locs) = true;

            [~,locs] = findpeaks(-AreaRateZvT(z,:));           
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
        % each loop over each one storing the rates (not propagation rate but area oscillation rate) and time points along
        % the labeled peak
        Nexp = max(Lexp,[],'all');
        ExpandMatIndividual = NaN(Nz,Nexp,2);
        for ii=1:Nexp
            bw = Lexp == ii;
            ind = find(bw);
            [row,col] = ind2sub(size(bw),ind);
            ratevec = AreaRateZvT(ind);
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
            ratevec = AreaRateZvT(ind);
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

            % subplot(2,2,1)
            % ADFe = AreaRateZvT;
            % ADFe(Lexp>0) = NaN;
            % imagesc(tvecMin,zvec*Zsf,ADFe)
            % title('Expansion Peaks')
            % 
            % subplot(2,2,2)
            % imagesc(tvecMin,zvec*Zsf,Lexp)
            % title('Expansion Peaks Labeled')
            % 
            % subplot(2,2,3)
            % DiffBasalCont = (ExpandMat(16,:,2)-ExpandMat(1,:,2));
            % histogram(DiffBasalCont)
            % title('Apical Expansion 4um delay time')

            subplot(1,3,1)
            imagesc(tvecMin,zvec*Zsf,AreaRateZvT);colorbar

            subplot(1,3,2)
            ADFc = AreaRateZvT;
            ADFc(Lcnt>0) = NaN;
            imagesc(tvecMin,zvec*Zsf,ADFc)
            colorbar;
            title('Contraction Peaks')
            
            
            subplot(1,3,3)
            imagesc(tvecMin,zvec*Zsf,Lcnt)
            title('Contraction Peaks Labeled')
%             
%             
%             subplot(2,3,6)
%             DiffBasalCont = (ContractMat(16,:,2)-ContractMat(1,:,2));
%             histogram(DiffBasalCont)
%             title('Apical Contraction 4um delay time')
%             
            pause;
        end



    end % cell
    
    
end % movie
%close all


figure (1)
%edges = -120:20:120;
%centers = -110:20:110;

centers = -(spf*12):spf:(spf*12);
d = diff(centers)/2;
edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];


Nedges = length(edges)-1;
ApicalContHist = NaN(22,Nedges);

posValuesA = [];
posValuesB = [];
negValuesA = [];
negValuesB = [];

for ii=1:22
    DiffBasalCont = spf*(ContractMatAll(ii+38,:,2)-ContractMatAll(ii,:,2));
    DiffBasalCont(isnan(DiffBasalCont)) = [];
    [N,~,bin]= histcounts(DiffBasalCont,edges,'Normalization', 'probability');

    ApicalContHist(ii,:) = N;

    % save positive and negative values for innerquartile range later
    if ii >= 1 && ii <=5
        positiveValues = DiffBasalCont(DiffBasalCont>0);
        posValuesA = horzcat(posValuesA,positiveValues);

        negativeValues = DiffBasalCont(DiffBasalCont<0);
        negValuesA = horzcat(negValuesA,negativeValues);

        positiveValues = [];
        negativeValues = [];
    elseif ii >= 18 && ii <= 22
        positiveValues = DiffBasalCont(DiffBasalCont>0);
        posValuesB = horzcat(posValuesB,positiveValues);

        negativeValues = DiffBasalCont(DiffBasalCont<0);
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
binsNegative = ApicalContHist(:,1:round(numColumns/2)-1);
binsPositive = ApicalContHist(:,round(numColumns/2)+1:end);
binsAll = horzcat(ApicalContHist(:,1:round(numColumns/2)-1),ApicalContHist(:,round(numColumns/2)+1:end));
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
yValsPosA = mean(ApicalContHist(apicalLayers,binPosA));

[iqrPosB,xValsPosB] = iqr(posValuesB);
binPosB = find(fix(centers)==fix(xValsPosB(1)) | fix(centers)==fix(xValsPosB(2)));
yValsPosB = mean(ApicalContHist(basalLayers,binPosB));

[iqrNegA,xValsNegA] = iqr(negValuesA);
binNegA = find(fix(centers)==fix(xValsNegA(1)) | fix(centers)==fix(xValsNegA(2)));
yValsNegA = mean(ApicalContHist(apicalLayers,binNegA));

[iqrNegB,xValsNegB] = iqr(negValuesB);
binNegB = find(fix(centers)==fix(xValsNegB(1)) | fix(centers)==fix(xValsNegB(2)));
yValsNegB = mean(ApicalContHist(basalLayers,binNegB));


% PLOTS
subplot(1,3,1)
imagesc(centers,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
title('\Delta z = 10 \mum')
xlabel('\Delta t (s)')
ylabel('Depth (\mum)')
colorbar
clim([0 0.12])
set(gca,'Fontsize',18)
xlim([-(spf*12) spf*12])
xticks([-80 -60 -40 -20 0 20 40 60 80])



subplot(1,3,2)
%plot(centers,mean(ApicalContHist(apicalLayers,:),1),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
a = area(centers(1:13),mean(ApicalContHist(apicalLayers,1:13),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%EdgeColor','none')
hold on
b = area(centers(13:25),mean(ApicalContHist(apicalLayers,13:25),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%'EdgeColor','none')
a.FaceAlpha = 0.5;
b.FaceAlpha = 0.5;
xlabel('\Delta t (s)')
ylabel('Probability')
set(gca,'Fontsize',18)
ylim([0 0.12])
xlim([-(spf*12) spf*12])
xticks([-80 -60 -40 -20 0 20 40 60 80])
grid on
%alpha(0.5)
c = area(centers(binNegA(1):binNegA(2)),mean(ApicalContHist(apicalLayers,binNegA(1):binNegA(2)),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%EdgeColor','none'))
d = area(centers(binPosA(1):binPosA(2)),mean(ApicalContHist(apicalLayers,binPosA(1):binPosA(2)),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%EdgeColor','none'))
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
% area(centers(1:13),mean(ApicalContHist(lateralLayers,1:13),1),'FaceColor',[0 0.52 0.54],'LineWidth',2)%EdgeColor','none')
% hold on
% area(centers(13:25),mean(ApicalContHist(lateralLayers,13:25),1),'FaceColor',[0.89 0.43 0],'LineWidth',2)%'EdgeColor','none')
% xlabel('\Delta t (s)')
% ylabel('Probability')
% set(gca,'Fontsize',18)
% ylim([0 0.12])
% xlim([-(spf*12) spf*12])
% xticks([-80 -60 -40 -20 0 20 40 60 80])
% grid on
% alpha(0.5)
% scatter(cogNegativeL,probNegL,75,'k','filled')
% scatter(cogPositiveL,probPosL,75,'k','filled')

subplot(1,3,3)
%plot(centers,mean(ApicalContHist(basalLayers,:),1),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
a = area(centers(1:13),mean(ApicalContHist(basalLayers,1:13),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%'EdgeColor','none')
hold on
b = area(centers(13:25),mean(ApicalContHist(basalLayers,13:25),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%'EdgeColor','none')
a.FaceAlpha = 0.5;
b.FaceAlpha = 0.5;
xlabel('\Delta t (s)')
ylabel('Probability')
set(gca,'Fontsize',18)
ylim([0 0.12])
xlim([-(spf*12) spf*12])
xticks([-80 -60 -40 -20 0 20 40 60 80])
grid on
%alpha(0.25)
c = area(centers(binNegB(1):binNegB(2)),mean(ApicalContHist(basalLayers,binNegB(1):binNegB(2)),1),'FaceColor',[0 0.52 0.54],'LineWidth',2);%EdgeColor','none'))
d = area(centers(binPosB(1):binPosB(2)),mean(ApicalContHist(basalLayers,binPosB(1):binPosB(2)),1),'FaceColor',[0.89 0.43 0],'LineWidth',2);%EdgeColor','none'))
%c.FaceAlpha = 0.75;
%d.FaceAlpha = 0.75;
scatter(cogNegativeB,probNegB,75,'k','filled')
scatter(cogPositiveB,probPosB,75,'k','filled')
%plot([xValsNegB(1),xValsNegB(1)], [0,yValsNegB(1)],'k','LineWidth',1)
%plot([xValsNegB(2),xValsNegB(2)], [0,yValsNegB(2)],'k','LineWidth',1)
%plot([xValsPosB(1),xValsPosB(1)], [0,yValsPosB(1)],'k','LineWidth',1)
%plot([xValsPosB(2),xValsPosB(2)], [0,yValsPosB(2)],'k','LineWidth',1)


% make boxplot of percentage positive vs negative
x1 = prcNegative;
x2 = prcPositive;
x3 = prcNegative(apicalLayers,:);
x4 = prcPositive(apicalLayers,:);
x5 = prcNegative(basalLayers,:);
x6 = prcPositive(basalLayers,:);
%x = [x1; x2; x3; x4; x5; x6];
x = [x1; x2];

g1 = repmat({'Negative'},numel(prcNegative),1);
g2 = repmat({'Positive'},numel(prcPositive),1);
g3 = repmat({'Negative Apical'},numel(prcNegative(apicalLayers,:)),1);
g4 = repmat({'Positive Apical'},numel(prcPositive(apicalLayers,:)),1);
g5 = repmat({'Negative Basal'},numel(prcNegative(basalLayers,:)),1);
g6 = repmat({'Positive Basal'},numel(prcPositive(basalLayers,:)),1);
%g = [g1; g2; g3; g4; g5; g6];
g = [g1; g2];

[~,p1] = ttest2(prcNegative,prcPositive)
[~,p1] = ttest2(prcNegative(apicalLayers,:),prcPositive(apicalLayers,:))
[~,p1] = ttest2(prcNegative(basalLayers,:),prcPositive(basalLayers,:))

[~,p1] = kstest2(prcNegative,prcPositive)
[~,p1] = kstest2(prcNegative(apicalLayers,:),prcPositive(apicalLayers,:))
[~,p1] = kstest2(prcNegative(basalLayers,:),prcPositive(basalLayers,:))
figure
boxplot(x,g,'Symbol','')
grid on
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontSize',20)
ylim([30 70])

end


% 
% subplot(1,4,1)
% %imagesc(edges(1:end-1),(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
% %imagesc(edges,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
% imagesc(centers,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
% title('\Delta z = 10 \mum','FontSize',16)
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Depth (\mum)','FontSize',16)
% colorbar
% %caxis([0 0.35])
% %xlim([-100 100])
% %xlim([-49 49])
% %xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
% %xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])
% 
% % subplot(1,4,2)
% % plot(edges(1:end-1),ApicalContHist([1 11 22],:))
% % legend('Apical 10','Lateral 10','Basal 10')
% % xlabel('\Delta t (s)','FontSize',16)
% % ylabel('Probability','FontSize',16)
% % title('+ = Apical leading')
% 
% 
% subplot(1,4,2)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(1,:),[0.1 -50 22 0.2 50 22]);
% %estimate = fitcurveMultiGauss1D_simple(centers,ApicalContHist(1,:),[0.1 -20 10 0.2 20 10]);
% plot(centers,ApicalContHist(1,:),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% xlim([-100 100])
% grid on
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% subplot(1,4,3)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(11,:),[0.1 -50 22 0.2 50 22]);
% %estimate = fitcurveMultiGauss1D_simple(centers,ApicalContHist(11,:),[0.1 -20 10 0.2 20 10]);
% plot(centers,ApicalContHist(11,:),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% xlim([-100 100])
% grid on
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% subplot(1,4,4)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(22,:),[0.1 -50 22 0.2 50 22]);
% %estimate = fitcurveMultiGauss1D_simple(centers,ApicalContHist(22,:),[0.1 -20 10 0.2 20 10]);
% plot(centers,ApicalContHist(22,:),'-o','Linewidth',2,'Color',[0.4940 0.1840 0.5560])
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% xlim([-100 100])
% grid on
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])

%end


        % set a threshhold on the minimun numbers of layers
%         minLayers= 10;
%         for ii=1:Nexp
%             sumL = sum(Lexp(:)==ii);
%             if sumL < minLayers
%                 Lexp(Lexp==ii) = 0;
%             end
%         end
%         for ii=1:Ncnt
%             sumL = sum(Lcnt(:)==ii);
%             if sumL < minLayers
%                 Lcnt(Lcnt==ii) = 0;
%             end
%         end
%         Nexp = max(Lexp,[],'all');
%         Ncnt = max(Lcnt,[],'all');
        

% function [] = fit_two_gaussians(x
% x = -100:100;
% y = [];
% 
% %G = @(L,A,R,W) (A.*exp(-(L-R).^2 / W.^2))
% F_error = @(L,A1,A2,R1,R2,W1,W2) sum(abs((A1.*exp(-(L-R1).^2 / W1.^2)+A2.*exp(-(L-R2).^2 / W2.^2) - y).^2))
% 
% ARW = fminsearch(@(z) F_error(x,z(1),z(2),z(3),z(4),z(5),z(6)),[1;2;3;]);
% 
% A=ARW(1);
% R=ARW(2);
% W=ARW(3);
%     
%     
% end

%         expandPeaksMat = zeros([size(ADF),4]);
%         cntrctPeaksMat = zeros([size(ADF),4]);
%         for z=1:Nz
%             [pks,locs,w,p] = findpeaks(ADF(z,:));   
%           
%             expandPeaksMat(z,locs,1) = pks;
%             expandPeaksMat(z,locs,2) = w;
%             expandPeaksMat(z,locs,3) = p;
%             expandPeaksMat(z,locs,4) = 1;
%             
%             [pks,locs,w,p] = findpeaks(-ADF(z,:));           
%             cntrctPeaksMat(z,locs,1) = pks;
%             cntrctPeaksMat(z,locs,2) = w;
%             cntrctPeaksMat(z,locs,3) = p;
%             cntrctPeaksMat(z,locs,4) = 1;
%             
%         end



% figure(1)
% fs = 16;
% subplot(2,2,1)
% idxFullDepth = all(isfinite(ExpandMatAll(:,:,2)));
% y = mean(ExpandMatAll(:,idxFullDepth,2),2,'omitnan');
% plot(zvec*Zsf,(y-y(1))*spf)
% xlabel('Depth (microns)','FontSize',fs)
% ylabel('Average Peak Timing (s)','FontSize',fs)
% title(['Expansion, N=',num2str(sum(idxFullDepth))],'FontSize',fs)
% grid on
% 
% subplot(2,2,2)
% idxFullDepth = all(isfinite(ContractMatAll(:,:,2)));
% y = mean(ContractMatAll(:,idxFullDepth,2),2,'omitnan');
% plot(zvec*Zsf,(y-y(1))*spf)
% xlabel('Depth (microns)','FontSize',fs)
% ylabel('Average Peak Timing (s)','FontSize',fs)
% title(['Contraction, N=',num2str(sum(idxFullDepth))],'FontSize',fs)
% grid on
% 
% subplot(2,2,3)
% y = mean(ExpandMatAll(:,:,1),2,'omitnan');
% plot(zvec*Zsf,y)
% xlabel('Depth (microns)','FontSize',fs)
% ylabel('Rate ','FontSize',fs)
% title('Expansion','FontSize',fs)
% grid on
% 
% subplot(2,2,4)
% y = mean(ContractMatAll(:,:,1),2,'omitnan');
% plot(zvec*Zsf,y)
% xlabel('Depth (microns)','FontSize',fs)
% ylabel('Rate ','FontSize',fs)
% title('Contraction','FontSize',fs)
% grid on
% 
% 
% figure(2)
% idx2 = all(isfinite(ContractMatAll(:,:,2)));
% A = ContractMatAll(:,idx2,2);
% [~,I] = sort(A);
% subplot(1,2,1)
% histogram(I(1,:))
% ylabel('Count')
% xlabel('Contraction Leading z-depth')
% idx2 = all(isfinite(ExpandMatAll(:,:,2)));
% A2=ExpandMatAll(:,idx2,2);
% [~,I] = sort(A2);
% subplot(1,2,2)
% histogram(I(1,:))
% ylabel('Count')
% xlabel('Expansion Leading z-depth')
% 
% edges = -60:10:60;
% figure(3)
% 
% subplot(2,2,1)
% DiffBasalCont = spf*(ContractMatAll(16,:,2)-ContractMatAll(1,:,2));
% histogram(DiffBasalCont,edges)
% title('Apical 4 um Contractions (positive=apical leading)','FontSize',fs)
% xlabel('Delay time to travel 4 microns (s)','FontSize',fs)
% ylabel('Count','FontSize',fs)
% 
% subplot(2,2,2)
% DiffBasalCont = spf*(ExpandMatAll(16,:,2)-ExpandMatAll(1,:,2));
% histogram(DiffBasalCont,edges)
% title('Apical 4 um Expansion (positive=apical leading)','FontSize',fs)
% xlabel('Delay time to travel 4 microns (s)','FontSize',fs)
% ylabel('Count','FontSize',fs)
% 
% subplot(2,2,3)
% DiffBasalCont = spf*(ContractMatAll(60,:,2)-ContractMatAll(45,:,2));
% histogram(DiffBasalCont,edges)
% title('Basal 4um Contractions  (positive=apical leading)','FontSize',fs)
% xlabel('Delay time to travel 4 microns (s)','FontSize',fs)
% ylabel('Count','FontSize',fs)
% 
% subplot(2,2,4)
% DiffBasalCont = spf*(ExpandMatAll(60,:,2)-ExpandMatAll(45,:,2));
% histogram(DiffBasalCont,edges)
% title('Basal 4 um Expansions (positive=apical leading)','FontSize',fs)
% xlabel('Delay time to travel 4 microns (s)','FontSize',fs)
% ylabel('Count','FontSize',fs)
% 
% 
% 
% figure(4)
% idxFullDepth = all(isfinite(ContractMatAll(:,:,2)));
% y = mean(ExpandMatAll(:,idxFullDepth,2),2,'omitnan');
% 
% [~,I] = min(ContractMatAll(:,:,2));
% idxBasalCleads = I>=55 & idxFullDepth;
% idxApicalCleads= I<= 5 & idxFullDepth;
% 
% n1 = sum(idxBasalCleads);
% n2 = sum(idxApicalCleads);
% 
% subplot(2,2,1)
% y = mean(ContractMatAll(:,idxBasalCleads,2),2,'omitnan');
% plot(zvec*Zsf,(y-y(end))*spf);hold on
% y2 = mean(ContractMatAll(:,idxApicalCleads,2),2,'omitnan');
% plot(zvec*Zsf,(y2-y2(1))*spf);hold off
% xlabel('Depth (microns)','FontSize',fs)
% ylabel('Average Peak Timing (s)','FontSize',fs)
% title(['Basal-leads N=',num2str(n1),', Apical-leads N=',num2str(n2)],'FontSize',fs)
% legend('B-leading Contraction','A-leading Contraction')
% grid on
% 
% 
% subplot(2,2,3)
% y = ContractMatAll(:,idxBasalCleads,2)*spf;
% y = y-repmat(y(end,:),60,1);
% plot(y);hold on
% 
% xlabel('Depth (layers)','FontSize',fs)
% ylabel('Average Peak Timing (s)','FontSize',fs)
% title('Basal Leading Contraction samples')
% 
% 
% grid on
% 
% 
% subplot(2,2,2)
% idxFullDepth = all(isfinite(ExpandMatAll(:,:,2)));
% [~,I] = min(ExpandMatAll(:,:,2));
% idxBasalEleads =  I>=55 & idxFullDepth;
% idxApicalEleads = I<= 5 & idxFullDepth;
% n1 = sum(idxBasalEleads);
% n2 = sum(idxApicalEleads);
% y = mean(ExpandMatAll(:,idxBasalEleads,2),2,'omitnan');
% plot(zvec*Zsf,(y-y(end))*spf);hold on
% y2 = mean(ExpandMatAll(:,idxApicalEleads,2),2,'omitnan');
% plot(zvec*Zsf,(y2-y2(1))*spf);hold off
% xlabel('Depth (microns)','FontSize',fs)
% ylabel('Average Peak Timing (s)','FontSize',fs)
% title(['Basal-leads N=',num2str(n1),', Apical-leads N=',num2str(n2)],'FontSize',fs)
% legend('B-leading Expansion','A-leading Expansion')
% grid on
% 
% 
% subplot(2,2,4)
% y = ExpandMatAll(:,idxBasalEleads,2)*spf;
% y = y-repmat(y(end,:),60,1);
% plot(y);hold on
% 
% xlabel('Depth (layers)','FontSize',fs)
% ylabel('Average Peak Timing (s)','FontSize',fs)
% title('Basal Leading Expansion samples')
% grid on
% 
% 
% 
% figure(5)
% edgesC = -15:0.1:15;
% Nc = NaN(numel(zvec*Zsf),numel(edgesC)-1);
% for ii=1:60;Nc(ii,:) = histcounts(ContractMatAll(ii,:,1),edgesC);end
% 
% subplot(2,2,1)
% imagesc(edgesC,zvec*Zsf,Nc)
% ylabel('Depth (layers)','FontSize',fs)
% xlabel('Contraction Rate (microns^2/min)','FontSize',fs)
% title('Contraction Rate Histogram')
% 
% edgesE = -15:0.1:15;
% Nc = NaN(numel(zvec*Zsf),numel(edgesC)-1);
% for ii=1:60;Ne(ii,:) = histcounts(ExpandMatAll(ii,:,1),edgesE);end
% 
% subplot(2,2,2)
% imagesc(edgesE,zvec*Zsf,Ne)
% ylabel('Depth (layers)','FontSize',fs)
% xlabel('Contraction Rate (microns^2/min)','FontSize',fs)
% title('Expansion Rate Histogram')
% 
% 
% diffTc = diff(ContractMatAll(:,:,2))*spf;
% sh = 15;
% diffTc = (ContractMatAll(1:(end-sh),:,2)-ContractMatAll((sh+1):end,:,2));
% edgesC = -10:1:10;
% Nc = NaN(size(diffTc,1),numel(edgesC)-1);
% for ii=1:size(diffTc,1);Nc(ii,:) = histcounts(diffTc(ii,:),edgesC,'Normalization', 'probability');end
% 
% subplot(2,2,3)
% imagesc(edgesC,zvec*Zsf,Nc)
% ylabel('Depth (layers)','FontSize',fs)
% xlabel('Contraction Delay over 4um (sec)','FontSize',fs)
% 
% 
% diffTe = (ExpandMatAll(1:(end-sh),:,2)-ExpandMatAll((sh+1):end,:,2));
% 
% edgesE = -10:1:10;
% Ne = NaN(size(diffTe,1),numel(edgesE)-1);
% for ii=1:size(diffTe,1);Ne(ii,:) = histcounts(diffTe(ii,:),edgesE,'Normalization', 'probability');end
% 
% subplot(2,2,4)
% imagesc(edgesE,zvec*Zsf,Ne)
% ylabel('Depth (layers)','FontSize',fs)
% xlabel('Expansion Delay over 4um (sec)','FontSize',fs)
% 
% 
% 
% figure(6)
% 
% cols = size(ExpandMatAll,2);
% StartStopE = NaN(2,cols);
% for ii=1:cols
%     vec = isfinite(ExpandMatAll(:,ii,1));
%     if any(vec)
%         Start = find(vec,1,'first');
%         Stop = find(vec,1,'last');
%         StartStopE(1,ii) = Start;
%         StartStopE(2,ii) = Stop;
%     end
% end
%     
% 
% cols = size(ContractMatAll,2);
% StartStopC = NaN(2,cols);
% for ii=1:cols
%     vec = isfinite(ContractMatAll(:,ii,1));
%     if any(vec)
%         Start = find(vec,1,'first');
%         Stop = find(vec,1,'last');
%         StartStopC(1,ii) = Start;
%         StartStopC(2,ii) = Stop;
%     end
% end
% 
% subplot(2,3,1)
% histogram(StartStopE(1,:)*Zsf)
% xlabel('Depth of Start (um)')
% title('Expansion Start Depth')
% 
% subplot(2,3,2)
% histogram(StartStopE(2,:)*Zsf)
% xlabel('Depth of Stop (um)')
% title('Expansion Stop Depth')
% 
% subplot(2,3,3)
% histogram(diff(StartStopE)*Zsf)
% xlabel('Total Depth')
% 
% 
% subplot(2,3,4)
% histogram(StartStopC(1,:)*Zsf)
% xlabel('Depth of Start (um)')
% title('Contraction Start Depth')
% 
% subplot(2,3,5)
% histogram(StartStopC(2,:)*Zsf)
% xlabel('Depth of Stop (um)')
% title('Contraction Stop Depth')
% 
% subplot(2,3,6)
% histogram(diff(StartStopC)*Zsf)
% xlabel('Total Depth')