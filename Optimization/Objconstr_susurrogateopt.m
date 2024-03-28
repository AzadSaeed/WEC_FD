function f = Objconstr_susurrogateopt(y,pdata)


obj      = OBJ_Calculations(y,pdata);

% Constraints: Only call from Objconstr_susurrogateopt if surrogate_opt is
% chosen as the algorith

if strcmpi(pdata.Opt.Solver,'surrogateopt')

    f.Fval   = obj;

    if strcmpi(pdata.General.CaseStudy,'PCL_opt') || strcmpi(pdata.General.CaseStudy,'L_opt') ...
            || strcmpi(pdata.General.CaseStudy,'LP_opt')|| strcmpi(pdata.General.CaseStudy,'P_opt')...
            || strcmpi(pdata.General.CaseStudy,'PC_opt')
        [c, ~]   = Constr(y,pdata);
        f.Ineq   = c;
    end

else

    f = obj;

end



end



