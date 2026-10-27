function [] = T2LeadingLayerAndJunctionalMyosin_LS(dataSets,intType)
close all

fs = 16;
zvec = 1:60;
maxZlayer = zvec(end);
zstepMic = 0.2635;

shiftSec = 30; 
maxlagSec = 60;


b2rcmap = b2r(-6,6);
b2rcmap = flip(b2rcmap,1);
lsf = 0.1;

SortIdxAll =[];
DeltaT = [];


Ntype = numel(dataSets);
for dat=1:Ntype
    data = dataSets{dat};
    Nmovies = numel(data);
    for mn=1:Nmovies
         spf = data(mn).SecPerFrame;
        shiftF = round(shiftSec/spf);
        maxlagF = round(maxlagSec/spf);
        maxlagZ = 40;
        lagvecZ = -maxlagZ:maxlagZ;
        lagvecF = -maxlagF:maxlagF;
        % load data
        cd(data(mn).Source)
        cd ..
        loaddata = load('gridAnalysis.mat');
        Tmatrix = loaddata.trackingMatrixZT(:,:,zvec);
        Lengths3D = Tmatrix(:,1:8:end,:)*lsf;
        Lengths3D(Lengths3D < 0.1) = NaN;
        [~,Nframes,Nlayers] = size(Lengths3D);
        
        
        if strcmpi('raw',intType)
            loadStruct = load('InterfaceBasedIntensityData.mat').InterfaceBasedIntensityData;
            InterfaceMyo = loadStruct.InterfaceIntensity_Myo(:,:,zvec);
        elseif strcmpi('tophat',intType)
            loadStruct = load('IntensityDataTOPHAT_ball8.mat').IntensityData;
            InterfaceMyo = loadStruct.InterfaceIntensity(:,:,zvec);
        else
            error('Unknown intType')
        end

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
            
            
            % get intensity data
            I_T1 = squeeze(InterfaceMyo(t1t3ints(int,1),:,:));
            I_T3 = squeeze(InterfaceMyo(t1t3ints(int,2),:,:));
            I = NaN(Nframes,Nlayers);
            finiteT1 = isfinite(I_T1);
            finiteT3 = isfinite(I_T3);
            I(finiteT1) = I_T1(finiteT1);
            I(finiteT3) = I_T3(finiteT3);
            IinterpF = I;
            Iinterp = inpaint_nans(I,2);
            Iinterp = Iinterp(tvec,:);
            IinterpF = filterImage3DpaddedEdges(Iinterp, 'Gauss', 3);
            
            finiteT3 = finiteT3(tvec,:);
            
            
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

            LinterpF = LinterpF';
            IinterpF = IinterpF';
            
            subplot(2,2,1)
            ax1 = gca;
            tvecMin = (tvec - tvec(1))*spf/60;
            imagesc(tvecMin,zvec*zstepMic,LinterpF);hold on
            %colormap(b2rcmap)
            caxis([-6 6])
            %plot(T2min,zvecMicrons,'LineWidth',2,'Color','g');
            M=contour(tvecMin,zvec*zstepMic,LinterpF,[0 0],'LineWidth',2,'Color','k');hold off
            colorbar
            ylabel('Z-depth (microns)','FontSize',fs)
            xlabel('Time (min)','FontSize',fs)
            set(gca,'FontSize',fs)
            colormap(ax1,b2rcmap)
            title('Length (microns)','FontSize',fs)
            
            subplot(2,2,2)
            ax2 = gca;

            imagesc(tvecMin,zvec*zstepMic,IinterpF);hold on
            M=contour(tvecMin,zvec*zstepMic,LinterpF,[0 0],'LineWidth',2,'Color','k');hold off
            colorbar
            ylabel('Z-depth (microns)','FontSize',fs)
            xlabel('Time (min)','FontSize',fs)
            set(gca,'FontSize',fs)
            colormap(ax2,'jet')
            caxis([min(IinterpF(:)) max(IinterpF(:))])
            title('Myosin Intensity','FontSize',fs)
            disp(num2str(int))
            
            
            LinterpF(~finiteT3') = NaN;
            IinterpF(~finiteT3') = NaN;
            subplot(2,2,3)
            [Lrate] = ShiftedGrad( LinterpF, shiftF);
            [Irate] = ShiftedGrad( IinterpF, shiftF);
            
            Lrate = normalize(Lrate,2);
            Irate = normalize(Irate,2);
            [ ccmat ] = crosscorrelation2D_varTemplateSize( Irate, Lrate, 20, maxlagF );
            [ acmat ] = crosscorrelation2D_varTemplateSize( Irate, Irate, 20, maxlagF );
            %c = xcorr2(Lrate,Irate);
            
            imagesc(lagvecF*spf,lagvecZ*zstepMic,ccmat);
            xlabel('\Delta T','FontSize',fs)
            ylabel('\Delta Z','FontSize',fs)
            title('Cross Correlation','FontSize',fs)
            colorbar
            %caxis([-0.5 0.5])
            grid on
            
            subplot(2,2,4)
            imagesc(lagvecF*spf,lagvecZ*zstepMic,acmat);
            xlabel('\Delta T','FontSize',fs)
            ylabel('\Delta Z','FontSize',fs)
            title('Myosin Rate Autco-correlation','FontSize',fs)
            colorbar
            %caxis([-0.5 0.5])
            grid on
            
            
            
           pause;

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


function [Signalsh] = ShiftedGrad( Signal, shift)



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

