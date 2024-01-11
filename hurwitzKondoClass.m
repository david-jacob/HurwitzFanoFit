
% hurwitzKondoClass
%
% 
% MatLAB script  for fitting Hurwitz-Fano / Frota-Fano lineshapes to experimental data
% (c) 2024 Robert Drost, Aalto University, Espoo, Finland

classdef hurwitzKondoClass

    %% Class properties:

    properties
        % Rquired input properties:
        xData = [];         % Typically bias voltage. Sets units of fit parameters.
        yData = [];         % Typically LIX data points 
        T = [];             % Temperature in K
        
        % Additional properties
        li = [];            % Lock-in modulation amplitude in units of xData

        % Optimisation parameters:
        a = 1;              % Scale factor
        Gamma = 1e-3;       % Encodes Kondo peak width
        phi = 0;            % Encodes Kondo peak asymmetry
        bgSlope = 0;        % Linear background slope
        bgIntercept = 0;    % Linear background y intercept
        biasOffset = 0;     % Bias offset

        % Optimisation result:
        fit = [];           % Fit curve

        % Optimisation options:
        useLi = true;       % Include lock-in broadening (true/false)
        model = 'hurwitz';  % Fitting model (hurwitz/frota)
        precision = 1e-6;   % Numerical precision for evaluating the Hurwitz Zeta function. Increase to improve speed.
        options = [];       % Optimisation options of lsqnonlin
        units = [];         % X-axis units (V/mV/uV)
        speed = true;       % Fast parameter initialisation from downsampled fit. Can improve speed considerably.



    % PROPERTIES end    
    end

    %% Object Functions

    methods
        
        % CONSTRUCTOR FUNCTION
        function obj = hurwitzKondoClass(xData, yData, temperature)
           
            % Validate inputs as numerical objects in the expected format:
            p = inputParser;
            vectorInput = @(x) isvector(x);
            scalarInput = @(x) isscalar(x);
            addRequired(p, 'xData', vectorInput);
            addRequired(p, 'yData', vectorInput);
            addRequired(p, 'temperature', scalarInput);
            parse(p, xData, yData, temperature);
            
            % Additional check for identical size of x and y data:
            if size(xData) == size(yData)
            
                % Set input properties
                obj.xData = xData;
                obj.yData = yData./max(yData);
                
                if sum(size(temperature)) == 2
                    obj.T  = temperature;
                else
                    obj.T = [];
                    fprintf('ERROR: Improperly formated input.\nTemperature must be a scalar.\nPlase check inputs.\n\n');
                end

                
                % Create optimisation options:
                obj.options = optimoptions('lsqnonlin');
                obj.options.Display = 'iter';
                
                % Guess xData input units:
                if max(abs(obj.xData)) > 0.1
                    obj.units = 'mV';
                    obj.Gamma = 1;
                else
                    obj.units = 'V';
                end
    
                fprintf(['Hurwitz-Kondo class initialised.\nUnits automatically set to ' obj.units '. Please adjust before fitting if necessary.\n']);

            else
                obj.xData = [];
                obj.yData = [];
                obj.T = [];
                fprintf('ERROR: Improperly formated input.\nxData and yData must be vectors of identical size.\nPlase check inputs.\n\n');
            end

            
        % CONSTRUCTOR FUNCTION end
        end



        % OBJECT FIT FUNCTION
        %  Fits the data stored in the object using the model specified by
        %  obj.model

        function obj = dofit(obj)

            fitStart = tic;

            if isempty(obj.T) || isempty(obj.xData) || isempty(obj.yData)
                fprintf('ERROR: Object xData, yData, and temperature must not be empty.\n');
                return;
            end
            
            % Define initial point
            xNaught = [obj.a, obj.Gamma/2.542, obj.phi, obj.biasOffset, obj.bgSlope, obj.bgIntercept];

            if strcmp(obj.model, 'hurwitz')
                fprintf('Fitting data with Hurwitz-Kondo model.\n');

                % Perform a quick fit on downsampled data for parameter
                % initialisation:

                if obj.speed
                    fprintf('Running fast parameter initialisation.\n')
                    
                    % Temporarily save xData an yData
                    storex = obj.xData;
                    storey = obj.yData;

                    % Downsample object xData and yData to 51 points:
                    obj.xData = linspace(min(obj.xData), max(obj.xData), 51);
                    obj.yData = interp1(storex, obj.yData, obj.xData);
                    
                    % Define optimisation function and perform fit:
                    fastOptiFunc = @(x) objectFitFunction(obj, x) - obj.yData;
                    obj.options.Display = 'off';
                    fitResult = lsqnonlin(fastOptiFunc, xNaught, [0, 0, 0, -1, -Inf, -Inf], [Inf, 10, 2*pi, 1, Inf, Inf], obj.options);

                    % Replace object parameters by fit output.
                    obj.a = fitResult(1);
                    obj.Gamma = fitResult(2)*2.542;
                    obj.phi = fitResult(3);
                    obj.biasOffset = fitResult(4);
                    obj.bgSlope = fitResult(5);
                    obj.bgIntercept = fitResult(6);
        
                    % Restore original xData and yData:
                    obj.xData = storex;
                    obj.yData = storey;
    
                    % Define new starting point for full optimisation:
                    xNaught = [obj.a, obj.Gamma/2.542, obj.phi, obj.biasOffset, obj.bgSlope, obj.bgIntercept];
                    obj.options.Display = 'iter';

                    fprintf('Downsampled fit converged. Proceeding to full fit.\n');
                end

            elseif strcmp(obj.model, 'frota')
                fprintf('Fitting data with Frota-Fano model.\n');
            end
            
            % Define optimisation problem
            optiFunc = @(x) objectFitFunction(obj, x) - obj.yData;

            % Perform fit
            fitResult = lsqnonlin(optiFunc, xNaught, [0, 0, 0, -1, -Inf, -Inf], [Inf, 10, 2*pi, 1, Inf, Inf], obj.options);

            fitEnd = toc(fitStart);

            % Replace object parameters by fit output.
            obj.a = fitResult(1);
            obj.Gamma = fitResult(2)*2.542;
            obj.phi = fitResult(3);
            obj.biasOffset = fitResult(4);
            obj.bgSlope = fitResult(5);
            obj.bgIntercept = fitResult(6);

            % Calculate optimisation function at fininal point:
            obj.fit = objectFitFunction(obj, [obj.a, obj.Gamma/2.542, obj.phi, obj.biasOffset, obj.bgSlope, obj.bgIntercept]);
            
            % Output warning in case no LI amplitude was specified, but LI
            % broadening should be included
            if obj.useLi && isempty(obj.li)
                fprintf('PLEASE NOTE: Lock-in broadening was requested, but no amplitude specified.\nCurrent fit result does NOT include lock-in broadening.\nPlease check inputs and run again, if necessary.\n\n')
            elseif obj.useLi && mean(diff(obj.xData)) > obj.li
                
            end

            % Plot object at final point if fit was successful
            if ~isempty(obj.fit)
                plot(obj);
            end

        % HURWITZ FIT FUNCTION end
        end


        % OBJECT PLOT FUNCTION:
        % Plots yData over xData with latest fit result:
        function obj = plot(obj)

            figure; hold;
            set(gca, 'XLimitMethod', 'tickaligned');
            set(gca, 'YLimitMethod', 'padded');
            plot(obj.xData, obj.yData);

            if ~isempty(obj.fit)
                plot(obj.xData, obj.fit, '--');
                text(0.1, 0.85, ['$\Gamma =$ ' num2str(obj.Gamma, '%.2e')], 'FontSize', 18, 'Units', 'normalized');
                text(0.1, 0.75, ['$\Phi =$ ' num2str(obj.phi, '%.2e')], 'FontSize', 18, 'Units', 'normalized');
                legend('Data', 'Fit', 'Location', 'ne');

            end

        % HURWITZ PLOT FUNCTION end
        end

    % METHODS end
    end

% CLASSDEF end
end






%% LOCAL FUNCTIONS

% objectFitFunction
% Constructs the specified optimisation function from object parameters
% Uses either the Hurwitz or Frota line shapes based on object properties
% Includes lock-in broadening on request
function x = objectFitFunction(obj, optiPars)
    
    % Define tau
    if strcmp(obj.units, 'V')
        tau = obj.T*pi*8.6173e-05;
    elseif strcmp(obj.units, 'mV')
        tau = obj.T*pi*8.6173e-02;
    elseif strcmp(obj.units, 'uV')
        tau = obj.T*pi*8.6173e1;
    else
        tau = obj.T*pi*8.6173e-05;
        fprintf('ATTENTION: Invalid x-axis units. Using V as default.');
    end

    % Calculate requested fitting model
    if strcmp(obj.model, 'hurwitz')
        x = optiPars(1).*sqrt(optiPars(2)/(8*tau)).*real(exp(1i*optiPars(3)).*hzeta(3/2, optiPars(2)/(2*tau) + 0.5 + 1i*(obj.xData - optiPars(4))./(2*tau), obj.precision)) + obj.xData.*optiPars(5) + optiPars(6);
        
    elseif strcmp(obj.model, 'frota')
        %x = optiPars(1).*imag(1i.*exp(1i*optiPars(3)).*sqrt(1i*optiPars(2)./(obj.xData + 1i*optiPars(2) - optiPars(4)))) + obj.xData.*optiPars(5) + optiPars(6);
        x = optiPars(1).*real(exp(1i*optiPars(3))./sqrt(1 + 1i*(obj.xData - optiPars(4))./optiPars(2))) + obj.xData.*optiPars(5) + optiPars(6);
        x = fermiDiracBroadening(obj.T, obj.xData, x);
    else
        fprintf('Invalid fitting model. Specify ''hurwitz'' or ''frota''');
    end

    % Include lock-in broadening if requested to do so and feasible
    if obj.useLi && ~isempty(obj.li)
        x = lockInBroadening(obj.li, obj.xData, x);
    end

end


% fermiDiracBroadening
% Fermi-Dirac broadening by convolution with resolution function
function y = fermiDiracBroadening(T, xData, yData)
    
    % Definitions:
    kB = 8.6173e-02;
    % Fermi-Dirac smearing function:
    beta = @(x) sech(x./(2.*kB.*T));
    % Compute convolution
    y = conv(yData, beta(xData.*1e3).^2./sum(beta(xData.*1e3).^2), 'same');
end


% Lock-in broadening by convolution with resolution function
function y = lockInBroadening(amplitude, xData, yData)
        
        % extend x-axis to suppress edge effects
        upx = linspace(min(xData) - 10*amplitude, max(xData) + 10*amplitude, ((max(xData) + 10*amplitude) - (min(xData) - 10*amplitude))/(0.1*amplitude));

        % Compute smearing function    
        liConv = real(sqrt(amplitude^2 - upx.^2));

        % Normalise area to unity to keep scale
        liConv = liConv./sum(liConv);

        % Interpolate yData on new x-axis
        upy = interp1(xData, yData, upx, 'linear', 'extrap');

        % Convolve with data
        y = conv(liConv, upy, 'same');

        % Interpolate on old x-axis and return
        y = interp1(upx, y, xData);
end

% Numerically compute Hurwitz Zeta function
function y = hzeta(s, a, precision)
    
    y = 0;
    k = 0;
    exitFlag = false;

    while ~exitFlag
        
        oldy = y;
        y = y + 1./(k + a).^s;
        k = k+1;

        if abs(oldy - y) < precision
            exitFlag = true;
        end

    end
end


        
