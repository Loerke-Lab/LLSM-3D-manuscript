function [YmatCellLayers] = FFTonCellAreaRates_CellStats_LSklar(dataSets)
% FFTonCellAreaRates_CellStats_LS  Performs the FFT (fast Fourier
% Transform) of the rate of cell area change. This function was used for
% the KLAR movies which had different parameters than the light sheet
% movies
%
% Input:       dataSets: structure that contains the following fields:                  
%                   -'ImageFileList': a list of file names for the
%                     cropped images.
%                   -'Source': the path to the cropped images.
%                   -'SecPerFrame': number of seconds per time frame
%
             
%
% Output:       YmatCelllayers:
%
%
% 08/14/22
fs = 16;

% rate time shift
shiftSec = 30;
% detrending order
DTorder = 1;
% flag for showing individual cell results (0 means doesn't show).

plotfft = 0;

% dute to different cell trajectories having different durations, we get
% different frequency vectors. So, we will interpolate the frequency
% vectors
interpf = 0.001:0.001:0.03;
Nf = length(interpf);
minDurationSec = 3*60;

% zvec is the depth vector to process
% zvec = 1:60;
% lsf = 0.1;
% zsf = 0.2635;
% Alayer = 1;
% Llayer = 32;
% Blayer = 60;

% the following lines are for Klar Y27
zvec =4:21;
lsf = 0.163;
zsf = 1;
Alayer = 1;
Llayer = 9;
Blayer = 17;

zvecFig = zvec-3;


YmatCellLayers = [];

Ntype = numel(dataSets);
for dat=1:Ntype
    
    data = dataSets{dat};
    Nmovies = numel(data);
    
    for mn=1:Nmovies
        
        spf = data(mn).SecPerFrame;
        shiftFrames = round(shiftSec/spf);
        minDurFrames = minDurationSec/spf;

        % move to the movie directory to load cell area data
        cd(data(mn).Source); cd ..
        loadGdata = load('GeometryData.mat');
        Area = loadGdata.Area(:,:,zvec)*lsf^2;
        Area = Area(:,1:99,:); %                    JUST FOR Y27 KLAR
        [Ncells,~,Nlayers] = size(Area);

        % get the rate of change of the area
        [AreaChange] = shiftedChange(Area,shiftFrames,spf);
        
        
        

        
        
        for ii=1:Ncells
            
            cellAreaChange = squeeze(AreaChange(ii,:,:))';
            
            [~,idx] = FiniteSequenceFromNans(mean(cellAreaChange),1);
            if numel(idx) < minDurFrames
                continue;
            end
            imagesc(cellAreaChange(:,idx));
            
            cellMat = NaN(Nlayers,Nf);
            for z=1:Nlayers
                AR = cellAreaChange(z,idx);
                
                
%                 [AR,~] = FiniteSequenceFromNans(AreaChange(ii,:,z),1);

                ARdt = detrend(AR,DTorder);
                [Y2,f,maxF] = analysisfft(ARdt,spf,plotfft);

                Y2i = interp1(f,Y2,interpf);

               



                cellMat(z,:) = Y2i;


            end % z
            
            YmatCellLayers = cat(3,YmatCellLayers,cellMat);

        end % cell
        
        
        

%         subplot(1,3,mn)
%         interpp = 1./interpf;
%         imagesc(interpf,zvec*0.5,YmatEachZ);
%         %xlabel('Period (s)','FontSize',16)
%         xlabel('Frequency (Hz)','FontSize',16)
%         ylabel('Depth (microns)','FontSize',16)
%         title('|Y(f)|','FontSize',16)
%         grid on


    end % movie
end % type
zvecMic = zvecFig*zsf;



figure(1)
subplot(1,2,1)
MovieMean = nanmean(YmatCellLayers,3);
imagesc(interpf,zvecMic,MovieMean)
colorbar
ylabel('Depth (microns)','FontSize',fs)
xlabel('Frequency (Hz)','FontSize',fs)
caxis([0.4 2.5])

Ncells = size(YmatCellLayers,3);



subplot(1,2,2)
x = [squeeze(YmatCellLayers(Alayer,10,:));squeeze(YmatCellLayers(Llayer,10,:));squeeze(YmatCellLayers(Blayer,10,:))];
g1 = repmat({'Apical'},Ncells,1);
g2 = repmat({'Lateral'},Ncells,1);
g3 = repmat({'Basal'},Ncells,1);
g = [g1; g2; g3];
boxplot(x,g,'Symbol','')
ylabel('FFT Amplitude','FontSize',fs)
title(['Ncells=',num2str(Ncells)],'FontSize',fs)
grid on
set(gca,'Ylim',[-1 5])

[~,p_AvL] = ttest2(squeeze(YmatCellLayers(Alayer,10,:)),squeeze(YmatCellLayers(Llayer,10,:)))
[~,p_AvB] = ttest2(squeeze(YmatCellLayers(Alayer,10,:)),squeeze(YmatCellLayers(Blayer,10,:)))
[~,p_LvB] = ttest2(squeeze(YmatCellLayers(Llayer,10,:)),squeeze(YmatCellLayers(Blayer,10,:)))



% figure(2)
% 
% 
% subplot(1,2,2)
% 
% MovieMean = nanmean(YmatCellLayers,3);
% imagesc(interpf,zvecMic,MovieMean)
% colorbar
% ylabel('Depth (microns)','FontSize',fs)
% xlabel('Frequency (Hz)','FontSize',fs)
% caxis([0.4 2.5])



 
end
% SUBFUNCTIONS



function [Y2,f,maxF] = analysisfft(signal,T,plotfft)
%analysisfft takes the fast fourier transform (fft) of the input signal
%with sampling period T, 


% The sample time of 0.5 seconds indicates, by the nyquist sampling
% theorem, that the highes frequency that can be measured will be 1 Hz.
% Thus for all interfaces we can measure frequencies between 0 and 1 Hz.
% Thus the independent variable, frequency, will run from 0 to 1 with a
% sampling that depends on the length of the signal.

%T = 0.5;                       % sample time (500 ms = 0.5 s)
Fs = 1/T;                       % sample rate
L = length(signal);             % Length of signal
t = (0:L-1)*T;                  % Time vector


NFFT = 2^nextpow2(L);           % Next power of 2 from length of signal
Y = fft(signal,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);



if nargin > 2 && plotfft
    % Plot single-sided amplitude spectrum.
    plot(f,2*abs(Y(1:NFFT/2+1))) 
    %title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('|Y(f)|','FontSize',16)
end

Y2= 2*abs(Y(1:NFFT/2+1));


[C,I] = max(Y2);
maxF = f(I);

end



function [SigChange] = shiftedChange(Signal,shift,tsf)

sh=shift;
tlen = size(Signal,2);

% cropped trajectory left
mat1 = Signal(:,1:(tlen-sh),:);
% cropped trajectory right
mat2 = Signal(:,(1+sh):tlen,:);
% change is difference between the two cropped trajectories
SigChange = (mat2-mat1)/(tsf*shift)*60;


end

