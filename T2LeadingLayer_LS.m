function [] = T2LeadingLayer_LS(dataSets)
close all

fs = 16;
zvec = 1:65;
maxZlayer = zvec(end);
zstepMic = 0.2635;
zvecMicrons = zvec*zstepMic;
Nz = numel(zvec);
b2rcmap = b2r(-6,6);
b2rcmap = flip(b2rcmap,1);
lsf = 0.1;

SortIdxAll =[];

DeltaTa = [];
DeltaT = [];


Ntype = numel(dataSets);
for dat=1:Ntype
    data = dataSets{dat};
    Nmovies = numel(data);
    for mn=1:Nmovies
        spf = data(mn).SecPerFrame;
        % load data
        cd(data(mn).Source)
        cd ..
        loaddata = load('gridAnalysis.mat');
        Tmatrix = loaddata.trackingMatrixZT(:,:,zvec);
        Lengths3D = Tmatrix(:,1:8:end,:)*lsf;
        Lengths3D(Lengths3D < 0.1) = NaN;
        [~,Nframes,Nlayers] = size(Lengths3D);

        % load type T1 interface numbers for T1 (1st column) and T3 (2nd column).
        loaddata = load('typeT1IntsV2.mat');
        t1t3ints = loaddata.typeT1Ints;
        Nints = size(t1t3ints,1);   
        for int=1:Nints

            L_T1 = squeeze(Lengths3D(t1t3ints(int,1),:,:));
            L_T3 = squeeze(Lengths3D(t1t3ints(int,2),:,:));
            L = NaN(Nframes,Nlayers);
            finiteT1 = isfinite(L_T1);
            finiteT3 = isfinite(L_T3);
            L(finiteT1) = L_T1(finiteT1);
            L(finiteT3) = -L_T3(finiteT3);


            % find starting lengths and ending lengths
            finite = isfinite(L);
            startL = find(any(finite,2),1,'first');
            endL   = find(any(finite,2),1,'last');
%             startFrame = max([startCells,startL]);
%             endFrame = min([endCells,endL]);
            tvec = startL:endL;

            % interpolate NaNs and if interpolated length is less than t.h. convert
            % to zeros
            Linterp =inpaint_nans(L,2); %methods 2, 4, and 5 worked best.
            Lnans = isnan(L);
            Lthresh = 5*lsf; % threshhold is 5 pixels
            T2nans = Lnans & (Linterp < Lthresh);
            Linterp(T2nans) = 0;
            Linterp = Linterp(tvec,:);

            LinterpF = filterImage3DpaddedEdges(Linterp, 'Gauss', 3);
            
%             [sx,~] = size(LinterpF);
%             
%             LinterpF= imresize(LinterpF,[sx,NlayersRS]);
            
            T2 = NaN(1,Nlayers);
            for j=1:Nlayers
                %lastpos = find(LinterpF(:,j)>0,1,'last');
                firstneg= find(LinterpF(:,j)<0,1,'first');
                if isempty(firstneg)
                    %lastpos = NaN;
                    firstneg = NaN;
                end
                T2(j) = firstneg;
            end
            
            if all(isnan(T2))
                continue;
            end
            
            includeExtrap = false;
            if any(isnan(T2))
                
                if includeExtrap
                    T2 = interp1(zvec(isfinite(T2)),T2(isfinite(T2)),zvec,'linear','extrap');
                    [~,SortIdx] = sort(T2);
                    SortIdxAll = [SortIdxAll;SortIdx];
                else
                    continue;
                end
            else
                [~,SortIdx] = sort(T2);
                SortIdxAll = [SortIdxAll;SortIdx];
            end
            
            
            T2min = T2*spf/60;
            
            maxT2 = max(T2min);
            [minT2,LL] = min(T2min);
            

            DeltaT = [DeltaT;maxT2 - minT2];

            
            
%             subplot(1,3,1)
%             tvecMin = (tvec - tvec(1))*spf/60;
%             imagesc(tvecMin,zvecMicrons,LinterpF');hold on
%             colormap(b2rcmap)
%             caxis([-6 6])
%             %plot(T2min,zvecMicrons,'LineWidth',2,'Color','g');
%             M=contour(tvecMin,zvecMicrons,LinterpF',[0 0],'LineWidth',2,'Color','k');hold off
%             colorbar
%             ylabel('Z-depth (microns)','FontSize',fs)
%             xlabel('Time (min)','FontSize',fs)
%             set(gca,'FontSize',fs)
%             
%             subplot(1,3,2)
%             plot(T2,'b')
%             
%             subplot(1,3,3)
%             plot(SortIdx,'b')
            
            
            
%            pause;

        end % T1 transition
    end % movie
end % type   

% nanrows = any(isnan(SortIdxAll),2);
% 
% SortIdxAll(nanrows,:) = [];
nrows = size(SortIdxAll,1);

% subplot(1,2,1)
% colors = jet(5);
% edges = [10 20 30 40 50]*0.5;
% for z = 1:NlayersRS
%     
%     sumZleads = sum(SortIdxAll==z);
%     ZleadsNorm = sumZleads/nrows;
%     
%     plot(edges,ZleadsNorm,'LineWidth',2,'Color',colors(z,:));hold on
%     
% end
% xlabel('Depth (microns)','FontSize',16)
% ylabel('Probability Leading','FontSize',16)
% legend('5','10','15','20','25')
fs =16;
%subplot(1,2,2)
figure(1)
bins = (1:10:maxZlayer)*zstepMic;
histogram(SortIdxAll(:,1)*zstepMic,bins,'Normalization','Probability')
xlabel('Depth (microns)','FontSize',16)
ylabel('Probability Leading','FontSize',16)
title(['N= ',num2str(nrows)],'FontSize',16)
%hc = colorbar;
% cb = linspace(1,Nz,16);
% set(hc, 'YTick',cb, 'YTickLabel',cb)




end % function

