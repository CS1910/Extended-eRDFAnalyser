function [mu,sigma,scaling,muE,sigmaE,scalingE,Plot] = ...
    Gaussfit_PDF2(x,Y,numPeaks,TagName)
% sigma is not the uncertainty on mu, but Gauss width

numCurves = length(Y);

rminphys = evalin('base', 'rminphys');

mu       = zeros(numCurves,numPeaks);
sigma    = zeros(numCurves,numPeaks);
scaling  = zeros(numCurves,numPeaks);
muE      = zeros(numCurves,numPeaks);
sigmaE   = zeros(numCurves,numPeaks);
scalingE = zeros(numCurves,numPeaks);

%% initialize plot
Col = [1,2,4,5,7,8,9,11,0,1,2,4,5,7,8,9,11,0];
Plot = figure('Name',['G(r) Gaussfits ', strrep(TagName, '$$^{\circ}$$', '°')], ...
    'units', 'normalized', ...
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
sgtitle(t, TagName, 'Interpreter', 'Latex');

yl = [1 1];
for iy = 1:numCurves
    
    y = Y{iy};
    
    %% find numPeaks peaks
    margin = 8;
    pind = [];
    firstpeak = 1;
    for i = (1+margin):(length(y)-margin)
        if length(pind) < numPeaks
            if y(i) > y(i-1) && y(i) > y(i+1)
                if y(i) > y(i-2) && y(i) > y(i+2)
                    if firstpeak == 1
                        if y(i) > 0
                            pind = [pind, i];
                            firstpeak = 0;
                        end
                    elseif firstpeak == 0
                        pind = [pind, i];
                        firstpeak = 0;
                    end
                end
            elseif firstpeak == 0
                if y(i) < y(i-1) && y(i) < y(i+1)
                    if y(i) < y(i-2) && y(i) < y(i+2)
                        pind = [pind, i];
                    end
                end
            end
        else
            break
        end
    end
    
    Delta = 0*(1:numPeaks);
    posPeak = (y(pind) > 0);
    i = 1;
    while i < length(posPeak)
        c = 0;  % number of same signs behind posPeak(i) with same sign as posPeak(i)
        while (i+c+1 <= length(posPeak)) && (posPeak(i) == posPeak(i+c+1))
            c = c + 1;
        end
        if c > 0
            for j = 0:c
                if (mod(j,2) == 1) && (i+j+1<=length(pind))
                    Delta(i+j) = (y(pind(i+j-1))+y(pind(i+j+1)))/2;
                elseif mod(j,2) == 1
                    Delta(i+j) = y(pind(i+j-1));
                end
            end
        end
        i = i + c + 1;
    end
    
    %% create fit data
    numFP = 4:margin;
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
    set(gca, 'Fontsize', 13.5);
    plot(x, y, 'Color', 'black', 'DisplayName', '$$G(r)$$')
    
    rc = evalin('base', 'rc');  % correction range
    
    title(['$$r_\mathrm{co}=\;$$',num2str(rc(iy)),'$$\,$$\AA'] ,'Interpreter', 'Latex');
    hold on
    for i = 1:length(pind)
        plot(xfit{i,end}(1:2:end), yfit{i,end}(1:2:end) ,'-o', ...
            'MarkerSize', 3, ...
            'Color', Color(Col(i)), ...
            'DisplayName', ['\textsf{Data + Fit ' int2str(i) '}'])
    end
    
    box on
    xlim([rminphys 10])
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
    
    for i = 1:length(pind)
        
        f = @(p,xdata)exp(-(xdata-p(1)).^2/(2*p(2)^2))*p(3)+Delta(i);
        
        p0 = [x(pind(i)), 0.2, 0.95*(y(pind(i)+3))-Delta(i)];
        
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
    
                
    for i = 1:length(pind)
        f = @(p,xdata)exp(-(xdata-p(1)).^2/(2*p(2)^2))*p(3)+Delta(i);
        
        xcurve = linspace(min(x),max(x),1000);
        p = [d1mu(i), d1sigma(i), d1scaling(i)];
        plot(xcurve, f(p,xcurve), ...
            'Color', Color( Col(i) ), ...
            'Handlevisibility', 'off');
            %'DisplayName', ['Fit ', int2str(i)]);
        
    end
    
    sh  = NaN*zeros(1,numPeaks-length(d1mu));  % spaceholder for filling the array
    mu(iy,:)       = [d1mu      , sh];
    muE(iy,:)      = [d1muE     , sh];
    sigma(iy,:)    = [d1sigma   , sh];
    sigmaE(iy,:)   = [d1sigmaE  , sh];
    scaling(iy,:)  = [d1scaling , sh];
    scalingE(iy,:) = [d1scalingE, sh];
    %
    %     mu(iy,:)       = d1mu(1:11)      ;
    %     muE(iy,:)      = d1muE(1:11)     ;
    %     sigma(iy,:)    = d1sigma(1:11)   ;
    %     sigmaE(iy,:)   = d1sigmaE(1:11)  ;
    %     scaling(iy,:)  = d1scaling(1:11) ;
    %     scalingE(iy,:) = d1scalingE(1:11);
    
    
end
% make last tile (doublicate of previous tile) to the legend only
for iy = 1:numCurves
    nexttile(iy);
    set(gca, 'Fontsize', 13);
    ylim([yl(1) yl(2)]);
end

nexttile(numCurves)
legend('NumColumns', 1, 'AutoUpdate', 'off', 'Location', 'EastOutside', ...
    'Interpreter', 'Latex');
%cla;
%set(gca, 'Visible', 'off');


end