function out = Calc_mean_cmt(Input, flag, SM)
% This function calculates the mean output from hydrodynamic surrogate
% models within the committee. To save computationl time and memory, mex
% codes are generated for the surrogate models. The ModelFlag "UseMex" is
% therefore necessary since "USEANN", which directly uses ANN models is not
% supported (the ANN models are not included in the repository).
% These models can be obtained by direclt contacting the code developers.
%
% Primary Controbuter: Saeed Azad, PHD


ModelFlag = 'USEMEX';

switch upper(ModelFlag)

    case 'USEMEX'

        CMT = ListMexFiles(flag,Input);
        out = mean(CMT,2);

    case 'USEANN'

        switch flag

            case 'A'
                CMT = SM.Net_1WEC{1,1};
            case 'B'
                CMT = SM.Net_1WEC{1,2};
            case 'Fe_r'
                CMT = SM.Net_1WEC{1,3};
            case 'Fe_i'
                CMT = SM.Net_1WEC{1,4};
            case 'A_11'
                CMT = SM.Net_2WEC{1,1};
            case 'A_12'
                CMT = SM.Net_2WEC{1,2};
            case 'B_11'
                CMT = SM.Net_2WEC{1,3};
            case 'B_12'
                CMT = SM.Net_2WEC{1,4};
            case 'Fe_r_11'
                CMT = SM.Net_2WEC{1,5};
            case 'Fe_i_11'
                CMT = SM.Net_2WEC{1,6};
        end

        % Extract the Model
        MDL  = CMT(1,:);

        % Use cell function to improve speed
        out_ = cellfun(@(c) c(Input),MDL,'UniformOutput',false);

        % Put in a matrix
        out_m = [out_{:}];

        % Find mean
        out = mean(out_m,2);

end


end



function CMT = ListMexFiles(flag, Input)

if strcmpi(flag,'A') || strcmpi(flag,'B') || strcmpi(flag,'Fe_r') || strcmpi(flag,'Fe_i') 

    for i = 1:10
        func_name = strcat('CMT_1WEC_QBC_',flag,'_', string(i),'_','CodeGen','(', 'Input',')');
        CMT(:,i) =eval(func_name);
    end
else
    for i = 1:5
        func_name = strcat('CMT_2WEC_QBC_',flag,'_', string(i),'_','CodeGen','(', 'Input',')');
        CMT(:,i) = eval(func_name);
    end
end


end

