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

%             subplot(2,2,1)
%             ADFe = AreaRateZvT;
%             ADFe(Lexp>0) = NaN;
%             imagesc(tvecMin,zvec*Zsf,ADFe)
%             title('Expansion Peaks')
%             
%             subplot(2,2,2)
%             imagesc(tvecMin,zvec*Zsf,Lexp)
%             title('Expansion Peaks Labeled')
            
%             subplot(2,2,3)
%             DiffBasalCont = (ExpandMat(16,:,2)-ExpandMat(1,:,2));
%             histogram(DiffBasalCont)
%             title('Apical Expansion 4um delay time')

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


% %%% Added 4/26/23 KB %%%
%centers = -60:10:60;
%centers = -120:20:120;
centers = -(spf*12):spf:(spf*12);
d = diff(centers)/2;
edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];
%%%

Nedges = length(edges)-1;
ApicalContHist = NaN(22,Nedges);

for ii=1:22
    DiffBasalCont = spf*(ContractMatAll(ii+38,:,2)-ContractMatAll(ii,:,2));
    DiffBasalCont(isnan(DiffBasalCont)) = [];
    [N,~]= histcounts(DiffBasalCont,edges,'Normalization', 'probability');
    
    ApicalContHist(ii,:) = N;
end


% %%% Added 5/1/23 KB %%%
% Nedges = length(edges)-1;
% ApicalContHist = NaN(7,Nedges);
% 
% for ii=1:7
%     DiffBasalCont = spf*(ContractMatAll(ii+53,:,2)-ContractMatAll(ii,:,2));
%     DiffBasalCont(isnan(DiffBasalCont)) = [];
%     [N,~]= histcounts(DiffBasalCont,edges,'Normalization', 'probability');
%     
%     ApicalContHist(ii,:) = N;
% end
% %%%

subplot(1,4,1)
%imagesc(edges(1:end-1),(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
%imagesc(edges,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
imagesc(centers,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
title('\Delta z = 10 \mum','FontSize',16)
xlabel('\Delta t (s)','FontSize',16)
ylabel('Depth (\mum)','FontSize',16)
colorbar
caxis([0 0.13])
%xlim([-49 49])
%xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])

% subplot(1,4,2)
% plot(edges(1:end-1),ApicalContHist([1 11 22],:))
% legend('Apical 10','Lateral 10','Basal 10')
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% title('+ = Apical leading')


subplot(1,4,2)
%estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(1,:),[0.1 -50 22 0.2 50 22]);
estimate = fitcurveMultiGauss1D_simple(centers,ApicalContHist(1,:),[0.1 -20 10 0.2 20 10]);
xlabel('\Delta t (s)','FontSize',16)
ylabel('Probability','FontSize',16)
%title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])

subplot(1,4,3)
%estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(11,:),[0.1 -50 22 0.2 50 22]);
estimate = fitcurveMultiGauss1D_simple(centers,ApicalContHist(11,:),[0.1 -20 10 0.2 20 10]);
xlabel('\Delta t (s)','FontSize',16)
ylabel('Probability','FontSize',16)
%title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])

subplot(1,4,4)
%estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(22,:),[0.1 -50 22 0.2 50 22]);
estimate = fitcurveMultiGauss1D_simple(centers,ApicalContHist(22,:),[0.1 -20 10 0.2 20 10]);
xlabel('\Delta t (s)','FontSize',16)
ylabel('Probability','FontSize',16)
%title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])



end


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