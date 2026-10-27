function [ExpandMatAll,ContractMatAll] = MyosinPropagationFindPeaksInt_LS(data,MovieNums,showEgs)
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
XYsf = 0.104;
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

    % load the angle shift
    loaddata = load('AngleShift.mat');
    AngleShift = loaddata.AngleShift;

    % load myosin intensity array (interior - M for medial) Myosin in channel 1
    loadStruct = load('InterfaceBasedIntensityData.mat').InterfaceBasedIntensityData;
    myoArray = loadStruct.InterfaceIntensity_Myo(:,:,zvec);

    % load tracking matrix with length and angles
    loadStruct = load('trackingMatrixZT_NodeFix.mat');
    Tmatrix = loadStruct.trackingMatrixZT_NodeFix(:,:,zvec);
    LengthMat = Tmatrix(:,1:8:end,:)*XYsf;
    AngleMat = Tmatrix(:,2:8:end,:)-AngleShift;
    idx = AngleMat < - 90;
    AngleMat(idx) = AngleMat(idx) + 180;
    idx = AngleMat > 90;
    AngleMat(idx) = AngleMat(idx) - 180;
    AngleMat = abs(AngleMat);

    Nints = size(myoArray,1); % number of cells to loop over

    for int = 1:Nints

        MyoZvT = squeeze(myoArray(int,:,:))'; % squeeze so we have depth (rows) versus time (columns)
        IntAng = squeeze(AngleMat(int,:,:))';
        IntLen = squeeze(LengthMat(int,:,:))';


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
        IntAng = IntAng(:,frameVec);
        IntLen = IntLen(:,frameVec);
        tvecMin = frameVec*spf/60;


        % median filtering is used to smooth over transient detection errors
        MyoZvTmed = MyoZvT;%medfilt2(AreaZvT,[5 5],'symmetric');
        
        % take the rate of change
        MyoZvT2 = MyoZvTmed(:,(DeltaTframes+1):end);
        MyoZvT1 = MyoZvTmed(:,1:(end-DeltaTframes));
        MyoRateZvT = (MyoZvT2 - MyoZvT1)/DeltaTminutes;

        % get the initial angle for the rate
        IntAngInitial = IntAng(:,1:(end-DeltaTframes));

        % get middle angle for the rate
        IntAng2 = IntAng(:,(DeltaTframes+1):end);
        IntAng1 = IntAng(:,1:(end-DeltaTframes));
        IntAngMiddle = (IntAng1+IntAng2)/2;
        
        % filter the rate so that findpeaks finds true peaks and not noise
        [MyoRateZvT] = filterImage3DpaddedEdges(MyoRateZvT, 'Gauss', [3,2,0]);

        % normalize myosin
        % normMethod = 1;
        % if normMethod == 1
        %     % normalize the entire matrix, 
        %     MyoRateZvT = (MyoRateZvT - mean(MyoRateZvT,'all'))/std(MyoRateZvT,[],'all');
        % 
        % end

        % store time-averaged junction angles to distinguish vertical from
        % horizontal interfaces
        % meanAngleVec = mean(IntAng,2,'omitnan'); 
        % meanAngleMat = cat(3,meanAngleMat,meanAngleVec);

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
        ExpandMatIndividual = NaN(Nz,Nexp,3);
        for ii=1:Nexp
            bw = Lexp == ii;
            ind = find(bw);
            [row,col] = ind2sub(size(bw),ind);
            ratevec = MyoRateZvT(ind);
            anglevec = IntAngInitial(ind);
            [B,I] = sort(row);
            if any(diff(B)==0) % branch exists, so skip to next
                continue;
            end
            RateVec = NaN(Nz,1);
            RateVec(B) = ratevec(I);
            TimeVec = NaN(Nz,1);
            TimeVec(B) = col(I);
            AngleVec = NaN(Nz,1);
            AngleVec(B) = anglevec(I);
            ExpandMatIndividual(:,ii,1) = RateVec;
            ExpandMatIndividual(:,ii,2) = TimeVec;
            ExpandMatIndividual(:,ii,3) = AngleVec;
        end


        Ncnt = max(Lcnt,[],'all');
        ContractMatIndividual = NaN(Nz,Ncnt,3);
        for ii=1:Ncnt
            bw = Lcnt == ii;
            ind = find(bw);
            [row,col] = ind2sub(size(bw),ind);
            ratevec = MyoRateZvT(ind);
            anglevec = IntAngInitial(ind);
            [B,I] = sort(row);
            if any(diff(B)==0) % branch exists, so skip to next
                continue;
            end
            RateVec = NaN(Nz,1);
            RateVec(B) = ratevec(I);
            TimeVec = NaN(Nz,1);
            TimeVec(B) = double(col(I));
            AngleVec = NaN(Nz,1);
            AngleVec(B) = anglevec(I);
            ContractMatIndividual(:,ii,1) = RateVec;
            ContractMatIndividual(:,ii,2) = TimeVec;
            ContractMatIndividual(:,ii,3) = AngleVec;
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


% get vertical interfaces and horizontal interfaces
expAngles = ExpandMatAll(:,:,3);
meanExpAngles = mean(expAngles,'omitnan');
idxExpVert = meanExpAngles < 20;
idxExpHorz = meanExpAngles > 70;
ExpandMatVert = ExpandMatAll(:,idxExpVert,:);
ExpandMatHorz = ExpandMatAll(:,idxExpHorz,:);

conAngles = ContractMatAll(:,:,3);
meanConAngles = mean(conAngles,'omitnan');
idxConVert = meanConAngles < 20;
idxConHorz = meanConAngles > 70;
ContractMatVert = ContractMatAll(:,idxConVert,:);
ContractMatHorz = ContractMatAll(:,idxConHorz,:);


%%% Added 4/26/23 KB %%%
% centers = -60:10:60;
% centers = -120:20:120;
centers = -(spf*8):spf:(spf*8);
d = diff(centers)/2;
edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];

% sizeHist = 22;
% sizeJump = 38; % 10 microns

% sizeHist = 41;
% sizeJump = 19; % 5 microns

% sizeHist = 30;
% sizeJump = 30; % 8 microns

% sizeHist = 37;
% sizeJump = 23; % 6 microns

sizeHist = 45;
sizeJump = 15; % 4 microns

%%% Added 5/1/23 KB %%%
Nedges = length(edges)-1;
ApicalExpVertHist = NaN(sizeHist,Nedges);
ApicalExpHorzHist = NaN(sizeHist,Nedges);
ApicalExpAllHist = NaN(sizeHist,Nedges);

for ii=1:sizeHist
    DiffBasalExpVert = spf*(ExpandMatVert(ii+sizeJump,:,2)-ExpandMatVert(ii,:,2)); % expansion (myosin increase)
    DiffBasalExpVert(isnan(DiffBasalExpVert)) = [];
    [N,~]= histcounts(DiffBasalExpVert,edges,'Normalization', 'probability');

    ApicalExpVertHist(ii,:) = N;

    DiffBasalExpHorz = spf*(ExpandMatHorz(ii+sizeJump,:,2)-ExpandMatHorz(ii,:,2)); % expansion (myosin increase)
    DiffBasalExpHorz(isnan(DiffBasalExpHorz)) = [];
    [N,~]= histcounts(DiffBasalExpHorz,edges,'Normalization', 'probability');

    ApicalExpHorzHist(ii,:) = N;

    DiffBasalExpAll = spf*(ExpandMatAll(ii+sizeJump,:,2)-ExpandMatAll(ii,:,2)); % expansion (myosin increase)
    DiffBasalExpAll(isnan(DiffBasalExpAll)) = [];
    [N,~]= histcounts(DiffBasalExpAll,edges,'Normalization', 'probability');

    ApicalExpAllHist(ii,:) = N;
end
%%%





figure
%subplot(1,4,1)
%imagesc(edges(1:end-1),(1:60)*0.2635,ApicalExpHist);title('Expansion from apical end')
%imagesc(edges,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
imagesc(centers,(1:60)*0.2635,ApicalExpVertHist);title('Contraction from apical end')
title('\Delta z = 4 \mum','FontSize',16)
xlabel('\Delta t (s)','FontSize',16)
ylabel('Depth (\mum)','FontSize',16)
colorbar
%caxis([0 0.15])
%xlim([-49 49])
%xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])
caxis([0 0.20]);
yticks([0 5 10 15])
set(gca,'FontSize',20)

figure
imagesc(centers,(1:60)*0.2635,ApicalExpHorzHist);title('Contraction from apical end')
title('\Delta z = 4 \mum','FontSize',16)
xlabel('\Delta t (s)','FontSize',16)
ylabel('Depth (\mum)','FontSize',16)
colorbar
%caxis([0 0.15])
%xlim([-49 49])
%xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])
caxis([0 0.20]);
yticks([0 5 10 15])
set(gca,'FontSize',20)

figure
imagesc(centers,(1:60)*0.2635,ApicalExpAllHist);title('Contraction from apical end')
title('\Delta z = 4 \mum','FontSize',16)
xlabel('\Delta t (s)','FontSize',16)
ylabel('Depth (\mum)','FontSize',16)
colorbar
%caxis([0 0.15])
%xlim([-49 49])
%xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])
caxis([0 0.20]);
yticks([0 5 10 15])
set(gca,'FontSize',20)

end

% subplot(1,4,2)
% plot(edges(1:end-1),ApicalExpHist([1 11 22],:))
% legend('Apical 10','Lateral 10','Basal 10')
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% title('+ = Apical leading')
% 
% subplot(1,4,2)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(1,:),[0.1 -50 22 0.2 50 22]);
% estimate = fitcurveMultiGauss1D_simple(centers,ApicalExpHist(1,:),[0.1 -20 10 0.2 20 10]);
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% subplot(1,4,3)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(11,:),[0.1 -50 22 0.2 50 22]);
% estimate = fitcurveMultiGauss1D_simple(centers,ApicalExpHist(11,:),[0.1 -20 10 0.2 20 10]);
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 
% subplot(1,4,4)
% %estimate = fitcurveMultiGauss1D_simple(edges(1:end-1),ApicalContHist(22,:),[0.1 -50 22 0.2 50 22]);
% estimate = fitcurveMultiGauss1D_simple(centers,ApicalExpHist(22,:),[0.1 -20 10 0.2 20 10]);
% xlabel('\Delta t (s)','FontSize',16)
% ylabel('Probability','FontSize',16)
% %title(['\mu_1= ',num2str(estimate(2)),', \mu_2= ',num2str(estimate(5)),', A_1/A_2= ',num2str(estimate(1)*estimate(3)/estimate(4)/estimate(6))])
% 

















