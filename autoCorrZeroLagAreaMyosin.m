function [] = autoCorrZeroLagAreaMyosin(data,MovieNum)

% this function plots the value at zero lag from the autocorrelations of the top layer area with every other layer
% also plots the value at zero lag from the autocorrelations of the top layer of myosin with every other layer


% INPUT:    


% OUTPUT:   just produces a plot with values of autocorrelations at zero lag across all depths
%           x-axis is depth, y-axis is autocorrelation value


% run the functions to produce autocorrelation matrices for area and myosin
[areaAutoCorr,maxLag] = cellAreaZLayers_LS(data,MovieNum,'AtoB');
[myosinAutoCorr,maxLag] = myosinZLayers_LS(data,MovieNum,'AtoB');

% the depth spacing and matrix of depth in microns
Zsf = 0.2635;
%Zsf = 1;
Nlayers = size(areaAutoCorr,3);
Zdepth = [1:Nlayers]*Zsf;

% find the autocorrelation values at zero lag for myosin and area across depths
zeroLag = maxLag+1;

% get the value at zero lag for both
zeroLagACArea = squeeze(areaAutoCorr(:,zeroLag,:))'; % rows are depth, columns are cells
zeroLagACMyo = squeeze(myosinAutoCorr(:,zeroLag,:))';

% take the average of all cells
avgACArea = mean(zeroLagACArea,2,'omitnan')';
avgACMyo = mean(zeroLagACMyo,2,'omitnan')';

% take the standard deviation of all cells
stdACArea = std(zeroLagACArea,[],2,'omitnan')';
stdACMyo = std(zeroLagACMyo,[],2,'omitnan')';

% get the number of cells
nCells = size(areaAutoCorr,1);

% get the standard error (std divided by square root of the sample size)
seACArea = (stdACArea/sqrt(nCells))';
seACMyo = (stdACMyo/sqrt(nCells))';

% plot the values at zero lag
figure
plot(Zdepth,avgACArea,'LineWidth',2,'Color',[0 0.4470 0.7410])
hold on
plot(Zdepth,avgACMyo,'LineWidth',2,'Color',[0.6350 0.0780 0.1840])
legend('Area Auto Corr','Myosin Auto Corr')
ylabel('Autocorrelation at Zero Lag')
xlabel('Depth (microns)')
set(gca,'Fontsize',20)
grid on


% plot the values at zero lag with standard error envelope
figure
xA = Zdepth';
yA = avgACArea';
eA = seACArea;
loA = yA - eA;
hiA = yA + eA;
hpA = patch([xA; xA(end:-1:1); xA(1)], [loA; hiA(end:-1:1); loA(1)], [0 0.4470 0.7410]);
alpha(0.25)
hold on;
hlA = line(xA,yA,'LineWidth',2);
grid on
set(hpA, 'edgecolor', 'none');
set(hlA, 'color', [0 0.4470 0.7410]);
xlabel('Depth (microns)','FontSize',12)
ylabel('Autocorrelation Zero Lag','FontSize',12)

xM = Zdepth';
yM = avgACMyo';
eM = seACMyo;
loM = yM - eM;
hiM = yM + eM;
hpM = patch([xM; xM(end:-1:1); xM(1)], [loM; hiM(end:-1:1); loM(1)], [0.6350 0.0780 0.1840]);
alpha(0.25)
hold on;
hlM = line(xM,yM,'LineWidth',2);
grid on
set(hpM, 'edgecolor', 'none');
set(hlM, 'color', [0.6350 0.0780 0.1840]);
xlabel('Depth (microns)','FontSize',12)
ylabel('Autocorrelation Zero Lag','FontSize',12)
set(gca,'Fontsize',20)
legend('','Area','','Myosin')

end

