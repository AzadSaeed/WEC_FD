function pdata = Wave_options(pdata, varargin)
% This function constructs the wave climate and modeling
% procedures. In addition to pdata, this function can take the region of
% the farm as an input in varargin. 
% 
% 
% Regular waves are defined based on a default significant wave height and
% wave period value.
% 
% Irregular waves are defined using Jonswap spectrum. Gauss quadrature is
% used for irregular waves. The irregular waves definition requires the
% availability of data for the specific site location. Currently, four
% locations are available:
%
% Alaska Coasts site ID NAWC24
% East Coast site ID NAEC8
% West Coast site ID N46229
% Pacific Islands site ID PI14
%
%
% Primary Contributor: Saeed Azad, PhD




switch upper(pdata.General.WaveType)

    case 'REGULAR'

        % Prescribed ragular wave significant wave height
        pdata.WaveData.Hs_regular   = 3;                                 

        % prescribed regular peack period
        pdata.WaveData.Te_regular   = 8.5;

        % Pass wave type to WaveData field
        pdata.WaveData.Type         = pdata.General.WaveType;

        % Create wave signal and power once outisde of the loop
        pdata                       = CallWaveSpectrum(pdata);

    case 'IRREGULAR'

        % [Hs Te] number of Gauss quadrature points in wave wave height
        % and Peak period dimensions to construct PDF
        n_GQ_SM             = [25, 20];
        pdata.WaveData.n_GQ = n_GQ_SM;


        if nargin == 1

            % Default is Alaska Coast
            pdata.WaveData.Region       = 'AlaskaCoasts';

        elseif nargin ==2

            % User-defined region
            pdata.WaveData.Region       = varargin{1,1};

        end


        
        switch pdata.WaveData.Region

            case 'AlaskaCoasts'

                % Alaskan Coasts Locations
                pdata.WaveData.site = 'NAWC24';

            case 'EastCoast'

                % East Coasts Locations
                pdata.WaveData.site = 'NAEC8';

            case 'PacificIslands'

                % Pacific Islands Locations
                pdata.WaveData.site = 'PI14';

            case 'WestCoast'

                % West Coast Locations
                pdata.WaveData.site = 'N46229';

        end

        % Climate model:
        % [] for historical,
        % 45 for mid century
        % 85 for end of century
        pdata.WaveData.climate = [];

        % Wave Data Period: historical or mid or old_data
        pdata.WaveData.WDP = 'historical';

        % Data for irregular wave form
        pdata                       = WaveInfo_func(pdata);
        
        % Number of regular waves to be used with superposition
        pdata.WaveData.nr           = 200;     

        pdata.WaveData.Type         = pdata.General.WaveType;
        
        % Jonswap
        pdata.WaveData.SpectrumType = 'JS'; 
        pdata.WaveData.BEMCount     = 200;

        % Create wave spectrum once outisde of the loop
        pdata                       = CallWaveSpectrum(pdata);
end

end