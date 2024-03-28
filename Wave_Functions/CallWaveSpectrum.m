function pdata = CallWaveSpectrum(pdata)


switch upper(pdata.WaveData.Type)

    case 'REGULAR'

        Hs             = pdata.WaveData.Hs_regular;
        Te             = pdata.WaveData.Te_regular;
        Omega          = (2*pi)/Te;
        S              = WaveSpectrum('REGULAR', Hs, Te, Omega, pdata);
        pdata.WaveData.S  = S;



    case 'IRREGULAR'

        switch pdata.WaveData.SpectrumType

            case 'JS'  

                % Local plot flag for demonstration
                pf         = 1;

                % Discretize angular frequency at nr points (for the superposition
                % of reglar waves)
                Omega          = (pdata.General.wMin:(pdata.General.wMax...
                    -pdata.General.wMin)/(pdata.WaveData.nr-1):pdata.General.wMax)';

                % Preallocate
                S   = cell(1,length(pdata.WaveData.nodes(:,1)));
                SP  = NaN(pdata.General.n_wec,length(Omega),length(pdata.WaveData.nodes(:,1)));
                DW  = NaN(pdata.General.n_wec,length(Omega),length(pdata.WaveData.nodes(:,1)));
                a_i = NaN(length(Omega),length(pdata.WaveData.nodes(:,1)));
                Ai  = NaN(pdata.General.n_wec,length(Omega),length(pdata.WaveData.nodes(:,1)));
                WavePower = NaN(1,length(pdata.WaveData.nodes(:,1)));
                for ii = 1:length(pdata.WaveData.nodes(:,1))

                    % Significant wave heights & wave periods
                    Hs         = pdata.WaveData.nodes(ii,1);
                    Te         = pdata.WaveData.nodes(ii,2);

                    % Perform Spertrum calculations
                    S{1,ii}    = WaveSpectrum('Jonswap', Hs, Te, Omega, pdata, pf);

                    % Extract spec and dw for power calculations
                    spec_        = S{1,ii}.spectrum.spectrum;
                    spec         = repmat(spec_,[1,pdata.General.n_wec])';
                    SP(:,:,ii)   = spec;
                    dw_          = S{1,ii}.spectrum.dOmega;
                    dw           = repmat(dw_,[1,pdata.General.n_wec])';
                    DW(:,:,ii)   = dw;

                    % Extract wave amplitude
                    a_i(:,ii)    = S{1,ii}.Regular_sup.Amplitude;
                    Ai(:,:,ii)   = repmat(a_i(:,ii),[1,pdata.General.n_wec])';
                   
                    % Power 
                    WavePower(1,ii) = S{1,ii}.spectrum.power;
                    
                end

        end

        pdata.WaveData.S  = S;
        pdata.WaveData.Ai = Ai;
        pdata.WaveData.DW = DW;
        pdata.WaveData.SP = SP;
        pdata.WaveData.WavePower = WavePower;
   

end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

