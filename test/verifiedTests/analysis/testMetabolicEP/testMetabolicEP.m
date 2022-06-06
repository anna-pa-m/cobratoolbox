% The COBRAToolbox: testMetabolicEP.m
%
% Purpose:
%     - run Expectation Propagation approximation to
%
% Authors:
%     - Anna Paola Muntoni

% load E. Coli core model
model = loadBiGGModel('e_coli_core');

% pre-process it (remove fixed fluxes, adjust lb and ub according to linear constraints)
pmodel = pre_processing(model);

% set the EP parameters
damping = 0.9;
precision = 1e-6;
maxit = 2000;
minvar = 1e-50;
maxvar = 1e50;
precision_lin = 1e-7;

% find Biomass flux index
idx_bm = find(matches(pmodel.rxns, 'BIOMASS_Ecoli_core_w_GAM'));
assert(~isempty(idx_bm))

% run Expectation Propagation for \beta \rightarrow +\infty
[muT0, sT0, aT0, dT0, avT0, vaT0, CovT0, t_EPT0] = MetabolicEPT0(full(pmodel.S), pmodel.b, ...
    pmodel.lb, pmodel.ub, damping, maxit, minvar, maxvar, precision, precision_lin, maxit);

% check average value for biomass flux
assert(abs(avT0(idx_bm) - 0.0418) < 1e-4)

% run Expectation Propagation for fixed \beta
warning('off','all')
Beta = 1e10;
[mu, s, a, d, av, va, Cov, t_EP]  = MetabolicEP(full(pmodel.S), pmodel.b, pmodel.lb, ...
    pmodel.ub,Beta, damping, maxit, minvar, maxvar, precision,  0, 0, 0, maxit);

% check average value for biomass flux
assert(abs(av(idx_bm) - 0.0432) < 1e-4)
fprintf('Test passed successfully\n')

