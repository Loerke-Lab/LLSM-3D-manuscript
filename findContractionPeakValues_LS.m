function [] = findContractionPeakValues_LS(data,MovieNums,showEgs)

% this function overlays contraction peaks on individual cell examples, then 
% calculates the value of the peak at each layer


% INPUTS:
%           data: structure that contains Source, ImageFileList, and
%                   SecPerFrame fields for each movie
%           MovieNums: vector of movie indices (of data) to include in
%                      analysis


%StartStop = [0 7];

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

% use the function b2r to generate a colormap, flip so it goes from red to blue
b2rcmap = flip(b2r(-6,6),1);

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
        %tvecMin = frameVec*spf/60;




        %%%%%
        % median filtering is used to smooth over transient detection errors
        AreaZvTmed = AreaZvT;%medfilt2(AreaZvT,[5 5],'symmetric');

        % take the rate of change
        AreaZvT2 = AreaZvTmed(:,(DeltaTframes+1):end);
        AreaZvT1 = AreaZvTmed(:,1:(end-DeltaTframes));
        AreaRateZvT = (AreaZvT2 - AreaZvT1)/DeltaTminutes;

        % filter the rate so that findpeaks finds true peaks and not noise
        [AreaRateZvT] = filterImage3DpaddedEdges(AreaRateZvT, 'Gauss', [3,2,0]);


        %%%%%
        % smooth out the area matrix
        % AreaAF = filterImage3DpaddedEdges(AreaZvT, 'Gauss', 2);
        % 
        % % calculate the rate of change and convert to microns^2/min
        % [AreaRateZvT] = RateOverDeltaT( AreaAF, DeltaTframes);
        % AreaRateZvT = AreaRateZvT/spf*60;



        tvecMin = (1:size(AreaRateZvT,2))*spf/60;
        

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
            subplot(1,3,1)
            imagesc(tvecMin,zvec*Zsf,AreaRateZvT);
            colormap(b2rcmap);
            colorbar
            maxL = max(abs(AreaRateZvT),[],'all');
            caxis(gca,[-maxL, maxL]);
    
            subplot(1,3,2)
            ADFc = AreaRateZvT;
            ADFc(Lcnt>0) = NaN;
            imagesc(tvecMin,zvec*Zsf,ADFc)
            colorbar;
            maxL = max(abs(AreaRateZvT),[],'all');
            caxis(gca,[-maxL, maxL]);
            title('Contraction Peaks')
            
            
            subplot(1,3,3)
            imagesc(tvecMin,zvec*Zsf,Lcnt)
            title('Contraction Peaks Labeled')
            colorbar
    
            %pause;
            %continue

            figure
            heatmap(tvecMin,zvec*Zsf,Lcnt)


            % make plot of contraction values at different depths
            contractPeak = 8;

            % get contraction values at that peak and z layer
            [rPeak,cPeak] = find(Lcnt==contractPeak);

            %[rPeakSort,idxSort] = sort(rPeak,'descend');
            [rPeakSort,idxSort] = sort(rPeak,'ascend');
            cPeakSort = cPeak(idxSort);

            peakRates = AreaRateZvT(Lcnt==contractPeak);
            peakRatesSort = peakRates(idxSort);
            zDepthSort = rPeakSort*Zsf;

            figure
            plot(zDepthSort,peakRatesSort,'-*r','LineWidth',2)
            %plot(zDepthSort(1:58),peakRatesSort(1:58),'-*r','LineWidth',2)
            %set(gca, 'XDir','reverse')
            grid on
            xlabel('Z Depth (microns)')
            ylabel('Peak Contraction Rate (microns/min)')
            set(gca,'Fontsize',20)
            ylim([-7 0])

            pause


        end


    end % cell
    
    
end % movie
%close all


%centers = -(spf*12):spf:(spf*12);
centers = -9.5:1:9.5; % makes an edge at zero
%centers = -10:1:10; % makes a center at zero
d = diff(centers)/2;
edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];

% Nedges = length(edges)-1;
% ApicalContHist = NaN(22,Nedges);
% 
% for ii=1:22
%     DiffBasalCont = (ContractMatAll(ii+38,:,1)-ContractMatAll(ii,:,1));
%     DiffBasalCont(isnan(DiffBasalCont)) = [];
%     [N,~]= histcounts(DiffBasalCont,edges,'Normalization', 'probability');
% 
%     ApicalContHist(ii,:) = N;
% end


% sizeHist = 22;
% sizeJump = 38; % 10 microns

% sizeHist = 30;
% sizeJump = 30; % 8 microns

sizeHist = 41;
sizeJump = 19; % 5 microns

%%%
Nedges = length(edges)-1;
ApicalContHist = NaN(sizeHist,Nedges);
cntNum = size(ContractMatAll,2);


for ii=1:sizeHist

    DiffBasalCont = NaN(1,cntNum);

    for cnt = 1:cntNum

        %if ContractMatAll(ii+sizeJump,cnt,2) >= ContractMatAll(ii,cnt,2)
        if ContractMatAll(ii+sizeJump,cnt,2) > ContractMatAll(ii,cnt,2)
            DiffBasalCont(:,cnt) = (ContractMatAll(ii+sizeJump,cnt,1)-ContractMatAll(ii,cnt,1));
        elseif ContractMatAll(ii+sizeJump,cnt,2) < ContractMatAll(ii,cnt,2)
            DiffBasalCont(:,cnt) = (ContractMatAll(ii,cnt,1)-ContractMatAll(ii+sizeJump,cnt,1));
        end

    end

    DiffBasalCont(isnan(DiffBasalCont)) = [];
    [N,~]= histcounts(DiffBasalCont,edges,'Normalization', 'probability');
    %[N,~]= histcounts(DiffBasalCont,edges);

    ApicalContHist(ii,:) = N;

    % prcPositive(ii,:) = (sum(DiffBasalCont>=0)/size(DiffBasalCont,2))*100;
    % prcNegative(ii,:) = (sum(DiffBasalCont<0)/size(DiffBasalCont,2))*100;


end
%%%


%subplot(1,4,1)
%imagesc(edges(1:end-1),(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
%imagesc(edges,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
imagesc(centers,(1:60)*0.2635,ApicalContHist);title('Contraction from apical end')
title('\Delta z = 10 \mum','FontSize',16)
xlabel('\Delta Rate','FontSize',16)
ylabel('Depth (\mum)','FontSize',16)
colorbar
caxis([0 0.12])
set(gca,'FontSize',20)
%xlim([-49 49])
%xticks([-54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2])
%xticks([-81.4 -67.8 -54.2 -40.7 -27.1 -13.6 0 13.6 27.1 40.7 54.2 67.8 81.4])

% get the percent of positive values vs the percent of negative values
numColumns = size(ApicalContHist,2);
%binsPositive = ApicalContHist(:,1:(numColumns/2));
%binsNegative = ApicalContHist(:,(numColumns/2)+1:end);

binsNegative = ApicalContHist(:,1:(numColumns/2));
binsPositive = ApicalContHist(:,(numColumns/2)+1:end);
prcPositive = (sum(binsPositive,'all')/sum(ApicalContHist,'all'))*100;
prcNegative = (sum(binsNegative,'all')/sum(ApicalContHist,'all'))*100;

% make boxplot of percentage of changes at each z depth
prcPositiveZ = (sum(binsPositive,2)./sum(ApicalContHist,2))*100;
prcNegativeZ = (sum(binsNegative,2)./sum(ApicalContHist,2))*100;

figure
x = [prcNegativeZ; prcPositiveZ];
g1 = repmat({'Negative Change'},sizeHist,1);
g2 = repmat({'Positive Change'},sizeHist,1);
g = [g1; g2];
boxplot(x,g,'Symbol','')
ylabel('Percent')
title('Contraction Rate Change \Delta z = 10 \mum')
grid on


x2 = reshape(x,sizeHist,[]);
hold on
plot(mean(x2),'dg')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontSize',20)
ylim([40 60])

[h,p] = kstest2(prcPositiveZ,prcNegativeZ)

end



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


    