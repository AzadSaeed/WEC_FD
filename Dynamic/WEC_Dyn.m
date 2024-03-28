function [P_total, pdata, varargout] = WEC_Dyn(pdata)

% This function calls the frequency-domain WEC dynamic models


% Calculate the output power [W]
[power_, pdata] = FD_Model(pdata);                       


switch upper(pdata.WaveData.Type)

    case 'REGULAR'

        % Estimated electrical power in W accounted for PCC efficiency
        Pe = power_.*pdata.General.PCC_eff;          

        % Report power as output since no probability is provided for regular waves [W]
        P_total = Pe;
        E       = Pe*24*365;                                                         
        p_mat   = [];
        Pae     = [];


    case 'IRREGULAR'

        p_mat      = NaN(pdata.WaveData.n_GQ(1,1),pdata.WaveData.n_GQ(1,2),pdata.General.n_wec);
        Pe         = NaN(pdata.WaveData.n_GQ(1,1),pdata.WaveData.n_GQ(1,2),pdata.General.n_wec);
        Pae        = NaN(pdata.WaveData.n_GQ(1,1),pdata.WaveData.n_GQ(1,2),pdata.WaveData.num_years,pdata.General.n_wec);
        P_total    = NaN(pdata.General.n_wec,pdata.WaveData.num_years);

        for c = 1:pdata.General.n_wec

            % Rearrange power into the appropriately-sized power matric [W]
            p_mat(:,:,c) = (reshape(wrev(power_(c,:)),[pdata.WaveData.n_GQ(1),pdata.WaveData.n_GQ(2)]));

            % Estimated electrical power in [W] accounted for PCC efficiency
            Pe(:,:,c) = p_mat(:,:,c).*pdata.General.PCC_eff;

            for  y = 1:pdata.WaveData.num_years
                Pae(:,:,y,c) = pdata.WaveData.Pr_mat(:,:,y).*Pe(:,:,c);
                P_total(c,y) = trapz(pdata.WaveData.nodes_Te,...
                    trapz(pdata.WaveData.nodes_Hs,squeeze(Pae(:,:,y,c)),1));
            end

        end

        % Total energy
        E = P_total*24*365;


        % Some plots
        if pdata.General.plotflag
            CreateSurfPlots(pdata, p_mat, Pae)
        end

end




if nargout > 2

    varargout = {power_, p_mat, Pe, Pae, E};

end


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function CreateSurfPlots(pdata,p_mat,Pae)

% define which year to plot
year = 1;

% define grid
Hs = pdata.WaveData.nodes_Hs;
Te = pdata.WaveData.nodes_Te;
[XX,YY] = ndgrid(Hs,Te);

commonFigureProperties;
hf = figure;
hf.Color = 'w';
surf(XX, YY, pdata.WaveData.Pr_mat(:,:,year))
xlabel('$H_s$')
ylabel('$T_e$')
zlabel('$Probability~matrix$')
savename  = strcat(pdata.General.plot_dir, filesep, 'Pr_matrix','.pdf');
exportgraphics(hf, savename)
close(hf)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
surf(XX, YY, p_mat(:,:,year))
xlabel('$H_s$')
ylabel('$T_e$')
zlabel('$Power~matrix$')
savename  = strcat(pdata.General.plot_dir, filesep, 'Power_matrix','.pdf');
exportgraphics(hf, savename)
close(hf)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
surf(XX, YY, Pae(:,:,year))
xlabel('$H_s$')
ylabel('$T_e$')
zlabel('$Power~matrix \times Probability~matrix$')
savename  = strcat(pdata.General.plot_dir, filesep, 'Probabilistic_Power_matrix','.pdf');
exportgraphics(hf, savename)
close(hf)

end
