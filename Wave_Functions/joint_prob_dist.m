function p_mat = joint_prob_dist(w_period,w_height,WaveData,pdata)


n_Gq  = WaveData.n_GQ(1,:);
p_mat = zeros(n_Gq(1),n_Gq(2),WaveData.num_years);

for y = 1:WaveData.num_years

    X              = [w_height{1,y}, w_period{1,y}];
    [p_mat_ln, xi]  = ksdensity(X,WaveData.nodes);
    p_mat_w        = p_mat_ln.*WaveData.Weights_mul;

    P_             = wrev(p_mat_w);
    p_mat_temp     = reshape(P_,[n_Gq(1),n_Gq(2)]);
    p_mat(:,:,y)   = p_mat_temp;

    hs_         = wrev(xi(:,1));
    Hs          = transpose(reshape(hs_,[n_Gq(1),n_Gq(2)]));

    te_         = wrev(xi(:,2));
    Te          =transpose(reshape(te_ ,[n_Gq(1),n_Gq(2)]));


    % Only plot for year 1
    if pdata.General.plotflag && y == 1

        commonFigureProperties;
        hf = figure;
        hf.Color = 'w';
        [XX,YY]     = ndgrid(transpose(Hs(1,:)),Te(:,1));
        contour(XX, YY, p_mat_temp,10)
        colormap(turbo)
        c = colorbar;
        clim([0 0.1])
        c.Label.String = 'Probability';
        view(0,90)
        xlabel('Hs')
        xlim([0,13])
        ylim([0, 25])
        ylabel('Te')
        zlabel('Probability density estimate')
        close all

    end

end
end







