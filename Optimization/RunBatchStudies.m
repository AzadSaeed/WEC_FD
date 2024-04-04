function pdata = RunBatchStudies(Regions,Studies,Hydro,IP, Ctrl, Location_Flag,...
    Solver, PowerLimit, WaveType, pdata)


for ii = 1:length(Regions)
    for jj = 1:length(Studies)
        for kk = 1:length(Hydro)
            for ll = 1:length(Ctrl)
                for mm = 1:length(IP)
                    for nn = 1:length(Location_Flag)
                        for oo = 1:length(Solver)
                            for pp = 1:length(PowerLimit)
                                for qq = 1:length(WaveType)

                                % Pass data to pdata structure
                                pdata.General.CaseStudy  = Studies{1,jj};

                                % Create directories to store the results
                                pdata = CreateDirectories(pdata);

                                % Precribe all problem data and problem settings
                                pdata = ProblemOptions(pdata,Regions{1,ii},...
                                    Hydro{1,kk},Ctrl{1,ll},IP{1,mm},...
                                    Location_Flag{1,nn}, Solver{1,oo},... ...
                                    PowerLimit{1,pp},WaveType{1,qq});

                                % Solve Problem
                                pdata = RunCases(pdata);

                                close all;

                                % To avoid having multiple interactive
                                % sessions
                                delete(gcp('nocreate'))

                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


end