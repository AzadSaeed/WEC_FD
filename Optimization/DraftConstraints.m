function c = DraftConstraints(D_min, D_max, Draft)

idx = 1;
c(idx,1) = D_min - Draft;

idx = idx+1;
c(idx,1) = Draft - D_max;

end