function [] = T2AlignedMeanLength_LS(dataSets,maxTimeMin)

od = cd;
fs = 16;
zvec = 1:60;
Nlayers = numel(zvec);
zstepMic = 0.2635;
zvecMicrons = zvec*zstepMic;
b2rcmap = b2r(-6,6);
b2rcmap = flip(b2rcmap,1);
close all
cmap = jet(Nlayers);
lsf = 0.1;

interpTimeVec = -maxTimeMin:0.2:maxTimeMin;
NframesInterp = length(interpTimeVec);
interpTimeMat = repmat(interpTimeVec,Nlayers,1);
Yq = repmat(zvec',1,NframesInterp);

DeltaTsec = 30;

LAligned3D = [];
RAligned3D = [];
T2all = [];

Ntype = numel(dataSets);
for dat=1:Ntype
    data = dataSets{dat};
    Nmovies = numel(data);
    for mn=1:Nmovies
        spf = data(mn).SecPerFrame;
        DeltaTf = round(DeltaTsec/spf);
        % load data
        cd(data(mn).Source)
        cd ..
        loaddata = load('gridAnalysis.mat');
        Tmatrix = loaddata.trackingMatrixZT(:,:,zvec);
        Lengths3D = Tmatrix(:,1:8:end,:)*lsf;
        Lengths3D(Lengths3D < 0.1) = NaN;
        Nframes = size(Lengths3D,2);

        % load type T1 interface numbers for T1 (1st column) and T3 (2nd column).
        loaddata = load('typeT1IntsV2.mat');
        t1t3ints = loaddata.typeT1Ints;
        Nints = size(t1t3ints,1);   
        for int=1:Nints

            L_T1 = squeeze(Lengths3D(t1t3ints(int,1),:,:))';
            L_T3 = squeeze(Lengths3D(t1t3ints(int,2),:,:))';
            L = NaN(Nlayers,Nframes);
            finiteT1 = isfinite(L_T1);
            finiteT3 = isfinite(L_T3);
            L(finiteT1) = L_T1(finiteT1);
            L(finiteT3) = -L_T3(finiteT3);


            % find starting lengths and ending lengths
            finite = isfinite(L);
            startL = find(any(finite,1),1,'first');
            endL   = find(any(finite,1),1,'last');
            tvec = startL:endL;
            NframesCrop = length(tvec);

            % interpolate NaNs and if interpolated length is less than t.h. convert
            % to zeros
            Linterp =inpaint_nans(L,2); %methods 2, 4, and 5 worked best.
%             Lnans = isnan(L);
%             Lthresh = 5*lsf; % threshhold is 5 pixels
%             T2nans = Lnans & (Linterp < Lthresh);
%             Linterp(T2nans) = 0;
            Linterp = Linterp(:,tvec);
            tvec = 1:numel(tvec);

            LinterpF = filterImage3DpaddedEdges(Linterp, 'Gauss', 2);
            
            medianT2 = find(median(LinterpF)<=0,1,'first');
            if isempty(medianT2)
                continue;
            end
           
            
            tvecNonAligend = tvec*spf/60;
            tvecAligned = (tvec - medianT2)*spf/60;
            
            if tvecAligned(1) > - maxTimeMin || tvecAligned(end) < maxTimeMin
                continue;
            end
            
            tvecAlignedMat = repmat(tvecAligned,Nlayers,1);
            Y = repmat(zvec',1,NframesCrop);
            
            LAligned = interp2(tvecAlignedMat,Y,LinterpF,interpTimeMat,Yq);
            
            anyTimeAllT1 = any(all(LAligned>0));
            anyTimeAllT3 = any(all(LAligned<0));
%             if ~anyTimeAllT1 || ~anyTimeAllT3
%                 continue;
%             end
            
            
            [Rates,centralTime] = centralDifference(-LinterpF,tvecAligned,DeltaTf,spf);
            
            centralTimeMat = repmat(centralTime,Nlayers,1);
            NframesCropRate = size(Rates,2);
            Y = repmat(zvec',1,NframesCropRate);
            
            RAligned = interp2(centralTimeMat,Y,Rates,interpTimeMat,Yq);
            
            
            T2 = NaN(1,Nlayers);
            for j=1:Nlayers
                %lastpos = find(LinterpF(:,j)>0,1,'last');
                firstneg= find(LAligned(j,:)<0,1,'first');
                if isempty(firstneg)
                    %lastpos = NaN;
                    firstneg = NaN;
                end
                T2(j) = firstneg;
            end
            
            if any(isnan(T2));continue;end
            [minval,ind] = min(T2);
            
            figure(1)
%             subplot(1,2,1)
%             imagesc(tvecNonAligend,zvecMicrons,LinterpF);
%             colormap(b2rcmap)
%             caxis([-6 6])
%             xlabel('Time (min)','FontSize',16)
%             ylabel('Depth (microns)','FontSize',16)
%             subplot(1,2,2)
            imagesc(interpTimeVec,zvecMicrons,LAligned);
            hold on;
            M=contour(interpTimeVec,zvecMicrons,LAligned,[0 0],'LineWidth',2,'Color','k');
            scatter(interpTimeVec(minval),ind*0.5,'k','filled');hold off
            colormap(b2rcmap);colorbar
            caxis([-6 6])
            xlabel('Time (min)','FontSize',16)
            ylabel('Depth (microns)','FontSize',16)
            title(['dataSets_mn',num2str(mn),', int',num2str(int)])
            
            pause;
            
            T2all = [T2all;ind*0.5];
            
            
            LAligned3D = cat(3,LAligned3D,LAligned);
            RAligned3D = cat(3,RAligned3D,RAligned);

            
%             T2min = T2*spf/60;
%             
%             tvecMin = (tvec - tvec(1))*spf/60;
%             imagesc(tvecMin,zvecMicrons,LinterpF');hold on
%             colormap(b2rcmap)
%             caxis([-6 6])
%             plot(T2min,zvecMicrons,'LineWidth',2,'Color','g');
%             M=contour(tvecMin,zvecMicrons,LinterpF',[0 0],'LineWidth',2,'Color','k');hold off
%             colorbar
%             ylabel('Z-depth (microns)','FontSize',fs)
%             xlabel('Time (min)','FontSize',fs)
%             set(gca,'FontSize',fs)
%             %pause;

        end % T1 transition
    end % movie
end % type   

meanL = mean(LAligned3D,3);
N = size(LAligned3D,3);

figure(1)

subplot(2,2,1)
imagesc(interpTimeVec,zvecMicrons,meanL);
ylabel('Z-depth (microns)','FontSize',fs)
xlabel('Time (min)','FontSize',fs)
colormap(b2rcmap)
colorbar
title(['Length, N=',num2str(N)])
set(gca,'FontSize',fs)

subplot(2,2,2)
 for j=1:Nlayers
    plot(interpTimeVec,meanL(j,:),'color',cmap(j,:),'LineWidth',2)
    hold on
end
hold off
grid on
ylabel('Length (microns)','FontSize',fs)
xlabel('Time (min)','FontSize',fs)
title(['Length, N=',num2str(N)])
set(gca,'FontSize',fs)
maxL = max(abs(meanL),[],'all');
set(gca,'Ylim',[-maxL maxL])



meanR = nanmean(RAligned3D,3);
N = size(RAligned3D,3);

subplot(2,2,3)
imagesc(interpTimeVec,zvecMicrons,meanR);
ylabel('Z-depth (microns)','FontSize',fs)
xlabel('Time (min)','FontSize',fs)
colormap(b2rcmap)
title(['Rate, N=',num2str(N)])
set(gca,'FontSize',16)

subplot(2,2,4)
 for j=1:Nlayers
    plot(interpTimeVec,meanR(j,:),'color',cmap(j,:),'LineWidth',2)
    hold on
end
hold off
grid on
ylabel('Rate (microns^2/min)','FontSize',fs)
xlabel('Time (min)','FontSize',fs)
title(['Rate, N=',num2str(N)])
set(gca,'FontSize',fs)
maxR = max(abs(meanR),[],'all');
set(gca,'Ylim',[-0 maxR])

cd(od)



figure(2);
LatT2 = LAligned3D(:,interpTimeVec==0,:);

Ledges = -4:1:4;
Lcenters = (Ledges(1:(end-1)) + Ledges(2:end))/2;
Nlbins = length(Lcenters);

LhistMat = NaN(Nlayers,Nlbins);
for z = zvec
    
    [N,~] = histcounts(squeeze(LatT2(z,:,:)),Ledges,'Normalization', 'probability');
    
    LhistMat(z,:) = N;
end
    

imagesc(Lcenters,zvecMicrons,LhistMat)
xlabel('Length at T=0','FontSize',fs)
ylabel('Depth (microns)','FontSize',fs)
title('Probability','FontSize',fs)
colorbar


figure(3)

bins = (1:5:60)*0.5;
histogram(T2all,bins)
xlabel('Depth (microns)','FontSize',16)
ylabel('Probability Leading','FontSize',16)
title(['N= ',length(T2all)],'FontSize',16)


end % function



function [Rates,centralTimeVec] = centralDifference(LengthVec,TimeVec,DeltaT,spf)



tlen = size(LengthVec,2);
    
idx1 = 1:(tlen-DeltaT);
idx2 = (1+DeltaT):tlen;
Rates = (LengthVec(:,idx2) - LengthVec(:,idx1))/DeltaT;
centralTimeVec = (TimeVec(:,idx1)+TimeVec(:,idx2))/2;

Rates = Rates/spf*60; % convert rates to min^{-1}

end

