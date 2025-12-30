function cfg = values()
%VALUES Parameter set for the Natasha configuration.

    cfg = struct();
    cfg.domainMinList = [0, 0];
    cfg.domainMaxList = [20, 20];
    cfg.domainGridSizeList = [401, 401];
    cfg.domainLabel = '\eta';

    p = struct();
    p.Pr = 30;
    p.lambda = -1.5;
    p.B = 0.5;
    p.E = 0.8;
    p.L1 = 0.5;
    p.L2 = 0.5;
    p.M = 0.2;
    p.S = 2.5;

    m = struct();
    m.phi_1 = 0.01; m.phi_2 = 0.01; m.phi_3 = 0.01;
    m.C_1 = 893;   m.rho_1 = 2720; m.k_1 = 222;  m.s_1 = 34.83e6;
    m.C_2 = 960;   m.rho_2 = 2810; m.k_2 = 173;  m.s_2 = 26.77e6;
    m.C_3 = 700;   m.rho_3 = 4907; m.k_3 = 3.7;  m.s_3 = 5.51e9;
    m.C_f = 1910;  m.rho_f = 884;  m.k_f = 0.144; m.s_f = 2.1e-12;

    m.s_nf   = m.s_f * ((m.s_3*(1+2*m.phi_3)+m.s_f*(1-2*m.phi_3))/(m.s_3*(1-m.phi_3)+m.s_f*(1+m.phi_3)));
    m.s_hnf  = m.s_nf * (m.s_2*(1+2*m.phi_2)+m.s_nf*(1-2*m.phi_2))/(m.s_2*(1-m.phi_2)+m.s_nf*(1+m.phi_2));
    m.s_thnf = m.s_hnf * (m.s_1*(1+2*m.phi_1)+m.s_hnf*(1-2*m.phi_1)*m.phi_1+m.s_2*m.phi_2)/(m.s_1*(1-m.phi_1)+m.s_hnf*(1+m.phi_1));
    m.k_nf   = m.k_f * ((2*m.k_f+m.k_3-2*m.phi_3*(m.k_f-m.k_3))/(m.k_3+2*m.k_f+m.phi_3*(m.k_f-m.k_3)));
    m.k_hnf  = m.k_nf*((2*m.k_nf+m.k_2-2*m.phi_2*(m.k_nf-m.k_2))/(m.k_2+2*m.k_nf+m.phi_2*(m.k_nf-m.k_2)));
    m.k_thnf = m.k_hnf*((m.k_1+2*m.k_hnf-2*m.phi_1*(m.k_hnf-m.k_1))/(m.k_1+2*m.k_hnf+m.phi_1*(m.k_hnf-m.k_1)));

    phiTotal = m.phi_1 + m.phi_2 + m.phi_3;

    n = struct();
    n.muRatio    = 1 / ((1 - phiTotal)^2.5);
    n.sigmaRatio = m.s_thnf / m.s_f;
    n.rhoRatio   = ((1-m.phi_2)*((1-m.phi_2)*((1-m.phi_3)*m.rho_f+m.phi_3*m.rho_3)+m.phi_2*m.rho_2)+m.phi_1*m.rho_1)/m.rho_f;
    n.kRatio     = m.k_thnf / m.k_f;
    n.rhoCpRatio = ((m.rho_1*m.C_1*m.phi_1)+(1-m.phi_1)*((1-m.phi_2)*((1-m.phi_3)*(m.rho_f*m.C_f))+m.phi_3*m.rho_2*m.C_2))/(m.rho_f*m.C_f);

    cfg.p = p;
    cfg.m = m;
    cfg.n = n;

    cfg.guesses = {
        @(x) [
            -exp(-x)
            exp(-x)
            exp(-x)
            exp(-x)
            exp(-x)
        ];
        @(x) natashaBranchTwoGuess(x, cfg.p)
    };
end

function guess = natashaBranchTwoGuess(x, ~)
    x = x(:).';
    onesRow = ones(1, numel(x));
    guess = [
        -exp(-x) + 0.9;
        -0.1 * onesRow;
        zeros(1, numel(x));
        zeros(1, numel(x));
        4 * onesRow
    ];
end
