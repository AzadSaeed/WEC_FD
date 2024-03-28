function pdata = dir_(pdata)

switch upper(pdata.General.CaseStudy)

    case 'PCL_OPT'

        switch upper(pdata.General.WaveType)

            case 'REGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_R;
                pdata.General.Solution_dir      = pdata.General.sol_dir_R;

            case 'IRREGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_IR;
                pdata.General.Solution_dir      = pdata.General.sol_dir_IR;
        end

    case 'C_OPT'

        switch upper(pdata.General.WaveType)

            case 'REGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_R;
                pdata.General.Solution_dir      = pdata.General.sol_dir_R;

            case 'IRREGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_IR;
                pdata.General.Solution_dir      = pdata.General.sol_dir_IR;
        end

    case 'L_OPT'

        switch upper(pdata.General.WaveType)

            case 'REGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_R;
                pdata.General.Solution_dir      = pdata.General.sol_dir_R;

            case 'IRREGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_IR;
                pdata.General.Solution_dir      = pdata.General.sol_dir_IR;
        end

    case 'P_OPT'

        switch upper(pdata.General.WaveType)

            case 'REGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_R;
                pdata.General.Solution_dir      = pdata.General.sol_dir_R;

            case 'IRREGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_IR;
                pdata.General.Solution_dir      = pdata.General.sol_dir_IR;
        end

    case 'LP_OPT'

        switch upper(pdata.General.WaveType)

            case 'REGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_R;
                pdata.General.Solution_dir      = pdata.General.sol_dir_R;

            case 'IRREGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_IR;
                pdata.General.Solution_dir      = pdata.General.sol_dir_IR;
        end

    case 'PC_OPT'

        switch upper(pdata.General.WaveType)

            case 'REGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_R;
                pdata.General.Solution_dir      = pdata.General.sol_dir_R;

            case 'IRREGULAR'

                pdata.General.plot_dir          = pdata.General.plot_dir_IR;
                pdata.General.Solution_dir      = pdata.General.sol_dir_IR;
        end

end


end