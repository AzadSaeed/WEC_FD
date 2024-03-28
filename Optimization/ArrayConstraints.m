function c = ArrayConstraints(xl,Radius,n_wec,SafeDist)

idx_  = 1:n_wec;
C     = nchoosek(idx_,uint16(2));
c     = zeros(length(C(:,1)),1);

for ii = 1: length(C(:,1))

    p_body  = C(ii,1);
    q_body  = C(ii,2);
    x_p     = xl(p_body,1);
    y_p     = xl(p_body,2);
    x_q     = xl(q_body,1);
    y_q     = xl(q_body,2);
    l_pq    = sqrt((x_p - x_q)^2 + (y_p - y_q)^2);
    c(ii,1) = 2*Radius + SafeDist - l_pq;


end

end