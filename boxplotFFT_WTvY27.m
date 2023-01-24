function [] = boxplotFFT_WTvY27(fftWT,fftY27)

YmatCellLayers = fftWT;
YmatCellLayers2 = fftY27;





Ncells1 = size(YmatCellLayers,3);
x1 = squeeze(YmatCellLayers(1,10,:));
x2 = squeeze(YmatCellLayers(30,10,:));
x3 = squeeze(YmatCellLayers(60,10,:));
g1 = repmat({'WTa'},Ncells1,1);
g2 = repmat({'WTl'},Ncells1,1);
g3 = repmat({'WTb'},Ncells1,1);
% x = [x1; x2; x3];
% g = [g1; g2; g3];


Ncells2 = size(YmatCellLayers2,3);
y1 = squeeze(YmatCellLayers2(1,10,:));
y2 = squeeze(YmatCellLayers2(30,10,:));
y3 = squeeze(YmatCellLayers2(60,10,:));
h1 = repmat({'Y27a'},Ncells2,1);
h2 = repmat({'Y27l'},Ncells2,1);
h3 = repmat({'Y27b'},Ncells2,1);
% y = [y1; y2; y3];
% h = [h1; h2; h3];


xy = [x1; y1; x2; y2; x3; y3];
gh = [g1; h1; g2; h2; g3; h3];

boxplot(xy,gh,'Symbol','')
ylabel('FFT Amplitude','FontSize',16)
%title(['Ncells=',num2str(Ncells)],'FontSize',fs)
grid on
set(gca,'Ylim',[-1 7])

[~,p_A] = kstest2(x1,y1)
[~,p_L] = kstest2(x2,y2)
[~,p_B] = kstest2(x3,y3)

end

