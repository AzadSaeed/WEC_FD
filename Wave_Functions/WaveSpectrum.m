function S = WaveSpectrum(SpectrumType, Hs, Tp, omega_, pdata,varargin)

% This function utilizes the used-defined inputs to characterize and 
% visualize regular or irregular wave.
% The analysis entails the computation of phase velocity, group velocit,
% wave power per unit crest (using exact equation and thus, avoiding
% deep/shallow water approximations.)
% For regular waves, the power is claculated as 
% p = (1/2)*rho*g*(amplitude^2)*cg
% The irregular wave is characterized through JONSWAP spectrum. 
% For irregular waves, the power is calculated as:
% p = rhp*g*integral(s(f)*cg*df)
% More detailed decsriptions are provided in the following link:
% https://wec-sim.github.io/WEC-Sim/master/theory/theory.html#irregular-waves
% Inputs
% SpectrumType ..... spectrum type, for now only Jonswap
% Hs           ............... significant Wave Height
% Tp           ............... peak period
% omega_       ............... radial frequency
% pdata        ............... porblem information
%
% Outputs
% S            ............... wave spectrum, including power 
%
% Contributor: Saeed Azad, PhD 

switch upper(SpectrumType)

    case 'REGULAR'

        theta = 0;

        % Nyquist frequency (highest observable frequency)
        % f_max = (L/2).*(Fs/L);          

        % Length of signal
        L = 6000;

        % x dimension 
        x = linspace(0,10,L);

        % Time vector
        t = linspace(0,70,L);

        % Water depth 
        [waterDepth, deepWater] = setWaterDepth(pdata.General.depth);

        % Calculate Wave Number for Larger Number of Frequencies
        wavenumber = calcWaveNumber(omega_,waterDepth,pdata.General.g,deepWater);                    
        
        % Significant Wave Height = sqrt(2)*Wave Height
        WaveHeight = (1/sqrt(2))*Hs;

        % Amplitude = Wave Height/2
        amplitude = WaveHeight/2;

        % Energy per unit horizontal sea surface area
        E = (1/8).*pdata.General.rho.*pdata.General.g.*(WaveHeight.^2);                     

        % Phase velocity 
        cp = CalculatePhaseVelocity(wavenumber,waterDepth,pdata.General.g);

        % Group Velocity [m/s]
        cg = CalculateGroupVelocity(wavenumber,waterDepth,pdata.General.g);

        % Calculate wave power per unit wave crest for regular waves - Do
        % not use estimations. Power calculated in [Watt/m] of wave crest
        power = (1/2)*pdata.General.rho*pdata.General.g*(amplitude^2)*cg;

        % Wave profile
        eta_f = @(t,x) amplitude*cos(wavenumber.*x - omega_.*(t) - theta);
        eta = eta_f(t,x);
    
        % Visualization
        if pdata.General.plotflag
            plotWaveProfile(t, eta, pdata)
        end

        S.power      = power; 
        S.E          = E;
        S.t          = t;
        S.eta        = eta; 
        S.amplitude  = amplitude;
        S.WaveHeight = WaveHeight;
        S.cg         = cg;
        S.cp         = cp;
        


    case 'JONSWAP'

        if nargin >=5
            % local plot flag 
            pf = varargin{1,1};   
        end

        % to be consistent with WECSIM
        BemCount = pdata.WaveData.BEMCount;

        % Angular frequency
        omega = (pdata.General.wMin:(pdata.General.wMax-...
            pdata.General.wMin)/(BemCount-1):pdata.General.wMax)'; 

        % Delta-omega
        dOmega = ones(BemCount,1).*(pdata.General.wMax-...
            pdata.General.wMin)./(BemCount-1);           

        % Wave Frequency
        frequency = omega/(2*pi);                                                                          

        % Define functions for coefficients and spectrum of pierson-Moskowitz  
        % B_ws coefficient for pierson-Moskowitz spectrum [1/(s^4)]
        bPM_fnc = @(Tp)(5/4)*(1/Tp)^(4);    

        % A_ws coefficient for pierson-Moskowitz spectrum [m^2/s^4]
        aPM_fnc = @(Hs,bPM) bPM*(Hs/2)^2;    

        % Wave Spectrum [m^2-s] for 'EqualEnergy' for Pierson-Moskowitx [m^2.s] or [m^2/Hz]
        fSpec_fnc = @(aPM,bPM,frequency)(aPM*frequency.^(-5).*exp(-bPM*frequency.^(-4)));                  

        % Evaluate coefficients and functions for pierson-Moskowitz
        % spectrum
        bPM = bPM_fnc(Tp);                                                                      
        aPM = aPM_fnc(Hs,bPM);   

        % [m^2.s]    
        fSpectrum = fSpec_fnc(aPM,bPM,frequency);                                                                                         
        [C,gammaAlpha] = JS_Spec_fnc(Hs,Tp,frequency);

        % Wave Spectrum [m^2.s] for Jonswap
        fSpectrum = C*fSpectrum.*gammaAlpha;

        % Wave Spectrum [m^2.s/rad]                                                               
        wspectrum = fSpectrum.*(2*pi); 
               
        % Water Depth
        [waterDepth, deepWater] = setWaterDepth(pdata.General.depth);

        % Calculate Wave Number for Larger Number of Frequencies [1/m]
        wavenumber = calcWaveNumber(omega,waterDepth,pdata.General.g,deepWater);                   
        
        % Group Velocity [m/s]
        cg = CalculateGroupVelocity(wavenumber,waterDepth,pdata.General.g);

        % Power per Unit Wave Crest: use full power equation 
        % Power in [Watt/m] of wave crest: [kg/m^3]*[m/s^2]*[m^2.s/rad]*[rad/s]*[m/s]
        power_vec = (pdata.General.rho*pdata.General.g.*wspectrum.*dOmega.*cg)./(2*pi)^2;
        power = sum(power_vec);

        % Create nr number of regular waves to use for superpositon        
        nr = pdata.WaveData.nr;                                                                    
        rng('default')
        theta = 2*pi.*rand(1,nr);
        
        % The following may be needed only if the Omega size is not consistent
        % dOmega_ = ones(nr,1).*(pdata.General.wMax-pdata.General.wMin)./(nr-1);
        % fr = omega_/(2*pi);
        % fSpectrum_      = fSpec_fnc(aPM,bPM,fr);
        % [C,gammaAlpha]  = JS_Spec_fnc(Hs,Tp,fr);
        % fSpectrum_      = C*fSpectrum_.*gammaAlpha;
        % wspectrum_      = fSpectrum_.*(2*pi);

        % Waveheight H = 2*A, A = sqrt(2*S*Delta_omega)(Dr. Herber's thesis - page 52)
        WaveHeight = 2.*sqrt(2*wspectrum.*(dOmega));                
        
        % Wave amplitude
        amplitude = WaveHeight/2;                                    

        % Potentially create the time signal from the spectrum
        TimeSignalFag = 0;

        % plot 
        % PlotJSSpectrum(omega, power_vec, pdata)
        
        if TimeSignalFag

            t = linspace(0,80,2000);
            x = linspace(0,10,2000);
            eta_ = zeros(nr,length(t));

            for i = 1:nr
                eta_(i,:) = (WaveHeight(i)./2).*cos(wavenumber(i).*x - ...
                    omega_(i).*t + theta(i));
            end

            E = (1/8).*pdata.General.rho.*pdata.General.g.*(WaveHeight.^2);
            cg = CalculateGroupVelocity(wavenumber,waterDepth,pdata.General.g);
            pwr = (pdata.General.rho*pdata.General.g.*wspectrum.*dOmega.*...
                cg)./(2*pi)^2;
            power_ = sum(pwr);
            eta = sum(eta_);

            if pf
                plotSpectrum(omega_,E,wspectrum,pdata);              
                plotWaveProfile(t, eta, pdata);
            end

        end


        
        % Power per Unit Wave Crest for the spectrum  
        % Power in [Watt/m] of wave crest
        S.spectrum.power_vec = power_vec;
        S.spectrum.power = power;                    

        % Spectrum in [m^2-s] 
        S.spectrum.fSpectrum = fSpectrum;                 

        % Spectrum in [m^2-s/rad]
        S.spectrum.spectrum = wspectrum;                 
        S.spectrum.dOmega = dOmega;

        % Outputs: Regular wave superposition information 
        % Number of regular waves for superposition
        S.Regular_sup.nr = nr;

        % WaveHeight
        S.Regular_sup.WaveHeight = WaveHeight;                       
        
        % Wave amplitudes
        S.Regular_sup.Amplitude = amplitude;

        if exist('eta','var')
            % Wave profile in time 
            S.Regular_sup.eta = eta;                      
        end


 end



end


%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

function [waterDepth, deepWater] = setWaterDepth(bemWaterDepth)

if strcmp(bemWaterDepth,'infinite')
    deepWater = 1;
    waterDepth = 200;
    fprintf('\tInfinite water depth specified in BEM and "waves.waterDepth" not specified in input file.\n')
    fprintf('Set water depth to 200m for visualization.\n')
else
    deepWater = 0;
    waterDepth = double(bemWaterDepth);
end
end

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

function [C gammaAlpha] = JS_Spec_fnc(Hs,Tp,frequency)
% This function calculates the required coefficient C and gammaAlpha in 
% order to construct the Jonswap spectrum from the availav=ble pierson 
% Moskowitz spectrum.

fp          = 1/Tp;                                                                                  % Peak frequency
siga        = 0.07;                                                                                  % Coefficient used in peak-shape parameter
sigb        = 0.09;                                                                                  % Coefficient used in peak-shape parameter
[lind,~]    = find(frequency<=fp);
[hind,~]    = find(frequency>fp);
gammaAlpha  = zeros(size(frequency));                                                                % nondimensional peak-shape parameter to the power of peak-shape parameter exponent
TpsqrtHs    = Tp/sqrt(Hs);

if TpsqrtHs <= 3.6
    gamma = 5;
elseif TpsqrtHs > 5
    gamma = 1;
else
    gamma = exp(5.75 - 1.15*TpsqrtHs);
end

gammaAlpha(lind) = gamma.^exp(-(frequency(lind)-fp).^2/(2*siga^2*fp^2));
gammaAlpha(hind) = gamma.^exp(-(frequency(hind)-fp).^2/(2*sigb^2*fp^2));
C = 1 - 0.287*log(gamma);

end

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

function cg = CalculateGroupVelocity(k,h,g)

kh  = k.*h;
cg = (1/2).*sqrt((g./k).*tanh(kh)).*(1+((2.*kh)./(sinh(2.*kh))));

end


%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

function cp = CalculatePhaseVelocity(k,h,g)

cp = sqrt((g/k)*tanh(k*h));

end

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

function plotSpectrum(omega_,E,wspectrum,pdata)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on

nr = pdata.WaveData.nr;
dOmega_ = ones(nr,1).*(pdata.General.wMax-pdata.General.wMin)./(nr-1);
E_scaled = (E/(pdata.General.rho*pdata.General.g))./dOmega_;

bar(omega_,E_scaled,'FaceColor',...
    [0.258823529411765	0.647058823529412	0.960784313725490],...
    'EdgeColor',[0.458823529411765	0.458823529411765	0.458823529411765])
plot(omega_,wspectrum,'LineWidth',linewidth, 'MarkerSize', 12,...
    'Color', C.red(10,:))
xlim([0 3.5]);
ax = gca;
ax.FontSize = fonttick;
xlabel('$\omega$', 'FontSize', fontlabel)
ylabel('$\textrm{Normalize~Spectrum}~[\frac{m^{2}s}{rad}]$', 'FontSize', fontlabel)

if strcmpi(pdata.General.WaveType,'Irregular')
    savename = strcat(pdata.General.plot_dir_IR, filesep,'Spectrum.pdf');
elseif strcmpi(pdata.General.WaveType,'Regular')
    savename = strcat(pdata.General.plot_dir_R, filesep,'Spectrum.pdf');
end
taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);

end

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
function plotWaveProfile(t, eta, pdata)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
plot(t,eta)          % time signal for the irregular wave
ax = gca;
ax.FontSize = fonttick;
xlabel('$t$', 'FontSize', fontlabel)
ylabel('$\textrm{Wave~profile}[m]$', 'FontSize', fontlabel)
savename = strcat(pdata.General.plot_dir, filesep,'IrregularWave.pdf');
taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);
close all

end
