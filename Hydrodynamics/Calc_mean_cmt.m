function out = Calc_mean_cmt(Input, flag, SM)

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