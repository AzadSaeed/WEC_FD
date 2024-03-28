function WaveData = gauss_q_pts_new(WaveData)


n_Gq_Hs = WaveData.n_GQ(1,1);
n_Gq_Te = WaveData.n_GQ(1,2);

a = WaveData.a;
b = WaveData.b;
c = WaveData.c;
d = WaveData.d;


% Define nodes and weights for guassian quadrature
[z1,w1] = lgwt(n_Gq_Hs,a,b);
[z2,w2] = lgwt(n_Gq_Te,c,d);



% Create node combinations without for loop
XZ       = {z1, z2};
C_1      = cell(size(XZ));
[C_1{:}] = ndgrid(XZ{:});
C_1      = cellfun(@(x) x(:),C_1,'un',0);
Cz       = [C_1{1,1:2}];   


% Create associated weight combinations without for loop
XW       = {w1, w2};
C_2      = cell(size(XW));
[C_2{:}] = ndgrid(XW{:});
C_2      = cellfun(@(x) x(:),C_2,'un',0);
CW_      = [C_2{1,1:2}];  
Cw       = CW_(:,1).*CW_(:,2);



WaveData.direction   = 0;       % Wave direction
WaveData.nodes       = Cz;      % Gauss quadrature node combination
WaveData.nodes_Hs    = z1;      % Gauss quadrature significant height nodes - vector
WaveData.nodes_Te    = z2;      % Gauss quadrature peak period nodes - vector   
WaveData.Weights_mul = Cw;      % Gauss quadrature weight combination
WaveData.Weights     = CW_; 


end

