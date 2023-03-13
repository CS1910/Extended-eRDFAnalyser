function [mu,sigma,scaling,muE,sigmaE,scalingE,Plot] = ...
    Gaussfit_Sq(x,Y,numPeaks,TagName)
% sigma is not the uncertainty on mu, but Gauss width

numCurves = length(Y);

mu       = zeros(numCurves,numPeaks);
sigma    = zeros(numCurves,numPeaks);
scaling  = zeros(numCurves,numPeaks);
muE      = zeros(numCurves,numPeaks);
sigmaE   = zeros(numCurves,numPeaks);
scalingE = zeros(numCurves,numPeaks);


%% initialize plot
Col = [1,2,4,5,7,8,9,11,0];
Plot = figure('Name',['S(q) Gaussfits ', strrep(TagName, '$$^{\circ}$$', '°')], ...
    'units', 'normalized',...
    'outerposition', [0.1 0 0.7 1], ...
    'NumberTitle', 'off');

factorials = factor(length(Y)+1);
if length(factorials) == 2
    t = tiledlayout(factorials(2),factorials(1));
else
    t = tiledlayout('flow');
end
%t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t, TagName, 'Interpreter', 'Latex');

yl = [1 1];
for iy = 1:numCurves
    
    y = Y{iy};
    
    %% find numPeaks peaks
    margin = 10;
    pind = [];
    for i = (1+margin):(length(y)-margin)
        if length(pind) < numPeaks
            if y(i) > y(i-1) && y(i) > y(i+1)
                if y(i) > y(i-2) && y(i) > y(i+2)
                    if y(i) > y(i-6) && y(i) > y(i+6)
                        pind = [pind, i];
                    end
                end
            end
        else
            break
        end
    end
    
    %% create fit data
    numFP = 6:margin;  % number of fitting points
    xfit = cell(numPeaks,length(numFP));
    yfit = cell(numPeaks,length(numFP));
    for j = 1:length(numFP)
        for i = 1:length(pind)
            xfit{i,j} = x((pind(i)-numFP(j)):(pind(i)+numFP(j)));
            yfit{i,j} = y((pind(i)-numFP(j)):(pind(i)+numFP(j)));
        end
    end
    
    %% subplot
    if iy == numCurves
        nexttile([1 2])
    else
        nexttile
    end
    set(gca, 'Fontsize', 13);
    plot(x, y+1, 'Color', 'black', 'DisplayName', '\textsf{$$S(q)$$}');
    
    rc = evalin('base', 'rc');  % correction range
    
    title(['$$r_\mathrm{co}=\;$$',num2str(rc(iy)),'$$\,$$\AA'], 'Interpreter', 'Latex');
    hold on
    for i = 1:length(pind)
        plot(xfit{i,end}(1:2:end), yfit{i,end}(1:2:end)+1, '-o', ...
            'MarkerSize', 3, ...
            'Color', Color(Col(i)), ...
            'DisplayName', ['\textsf{Data + Fit ' int2str(i) '}']);
    end
    
    box on
    xlim([1 10])
    % same y limits
    cylim = get(gca, 'ylim');
    if cylim(1) < yl(1)
        yl(1) = cylim(1);
    end
    if cylim(2) > yl(2)
        yl(2) = cylim(2);
    end
    
    
    %% fit % p1=mu, p2=sigma, p3=scaling
    %f = @(p,xdata)exp(-(xdata-p(1)).^2/(2*p(2)^2))*p(3)+Delta(i);
    
    d1mu       = [];
    d1muE      = [];
    d1sigma    = [];
    d1sigmaE   = [];
    d1scaling  = [];
    d1scalingE = [];
    
    % Gaussian with baseline at -0.05
    f = @(p,xdata)exp(-(xdata-p(1)).^2/(2*p(2)^2))*p(3)-0.05;
    
    for i = 1:length(pind)
        
        
        p0 = [x(pind(i)), 0.15, 0.95*y(pind(i)+3)+0.05];
        
        dmu       = [];
        dmuE      = [];
        dsigma    = [];
        dsigmaE   = [];
        dscaling  = [];
        dscalingE = [];
        for j = 1:length(numFP)
            [p,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(f,...
                p0, xfit{i,j}, yfit{i,j});
            ci = nlparci(p, resid, 'jacobian', J);
            dmu       = [dmu, p(1)];
            dmuE      = [dmuE, abs(p(1)-ci(1,1))];
            dsigma    = [dsigma, abs(p(2))];
            dsigmaE   = [dsigmaE, abs(p(2)-ci(2,1))];
            dscaling  = [dscaling, p(3)];
            dscalingE = [dscalingE, abs(p(3)-ci(3,1))];
        end
        d1mu       = [d1mu, mean(dmu)];
        d1muE      = [d1muE, std(dmu)];  % /sqrt(length(numFP))];
        d1sigma    = [d1sigma, mean(dsigma)];
        d1sigmaE   = [d1sigmaE, std(dsigma)];  % /sqrt(length(numFP))];
        d1scaling  = [d1scaling, mean(dscaling)];
        d1scalingE = [d1scalingE, std(dscaling)];  % /sqrt(length(numFP))];
        
    end
    
    % include in plot
    for i = 1:length(pind)
        xcurve = linspace(min(x),max(x),1000);
        p = [d1mu(i), d1sigma(i), d1scaling(i)];
        plot(xcurve, f(p,xcurve)+1, 'Color', Color(Col(i)), ...
            'HandleVisibility', 'off');
            %'DisplayName', ['Fit ', int2str(i)]);
    end
    
    % plot reference line at Sq=1
    xlims = get(gca, 'xlim');
    hold on
    plot([0 xlims(2)], [1 1], 'k', 'HandleVisibility', 'off');
    
    
    mu(iy,:)       = d1mu(:);
    muE(iy,:)      = d1muE(:);
    sigma(iy,:)    = d1sigma(:);
    sigmaE(iy,:)   = d1sigmaE(:);
    scaling(iy,:)  = d1scaling(:);
    scalingE(iy,:) = d1scalingE(:);
    
end
% make last tile (doublicate of previous tile) to the legend only
for iy = 1:numCurves
    nexttile(iy);
    set(gca, 'Fontsize', 11);
    ylim([yl(1) yl(2)]);
end

nexttile(numCurves)
legend('NumColumns', 1, 'AutoUpdate', 'off', 'Location', 'EastOutside', ...
    'Fontsize', 14, 'Interpreter', 'Latex');
% 
% % make last tile (doublicate of previous tile) to the legend only
% for iy = 1:numCurves+1
%     nexttile(iy);
%     set(gca, 'Fontsize', 11);
%     ylim([yl(1) yl(2)]);
% end
% nexttile(numCurves+1)
% set(gca, 'Fontsize', 14);
% legend('show', 'AutoUpdate', 'off');
% set(gca, 'Visible', 'off');
% cla;
end