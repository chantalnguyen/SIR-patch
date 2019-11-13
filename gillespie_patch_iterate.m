function [ t, x, rxnss, rxn_counts, XX ] = gillespie_patch_iterate( Npatch, s0, i0, n0, rep, gamma, l, thetastar, num_iter,patch_width)
%% Implementation of the Gillespie algorithm for N coupled population SIR model
%  using dynamically updated coupling matrix determined by von Mises
%  probability distribution with variance related to mean square
%  displacement
%   Usage:
%       [t, x, rxns, rxn_count] = gillespie_patch( tspan, Npatch, s0, i0, n0, rep, gamma, thetastar, num_iter )
%
%   Returns:
%       t:              cell of time vectors (Nreaction_events x num_iter)
%       x:              cell of averaged species amounts (Nreaction_events x Nspecies x num_iter)
%       rxnss:          reactions that occurred x num_iter
%       rxn_counts:     number of reactions undergone x num_iter
%       XX:             all species amounts (Nreaction_events x Nspecies x Npatch x num_iter)
%
%   Required:
%
%       Npatch:         Number of patches.
%
%       s0:             Initial number susceptible (Npatch x 1 vector).
%
%       i0:             Initial number infected (Npatch x 1 vector).
%
%       n0:             Total population of each patch (Npatch x 1 vector).
%
%       rep:            Basic reproductive number (beta/gamma)
%
%       gamma:          Recovery rate
%
%       l:              Step size
%
%       thetastar:      Maximum turning angle
%
%       num_iter:       Number of times to run simulation
%
%
%   Reference:
%       Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%       Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%
%       Nezar Abdennur, 2012 <nabdennur@gmail.com>
%       Dynamical Systems Biology Laboratory, University of Ottawa
%       www.sysbiolab.uottawa.ca
%
%   Created: 2012-01-19
%   Initially modified: 2019-10-23
%% Initialize
tspan = [0, 10000000];
if nargin < 10
	patch_width = 30; % width of each patch
end
[xx,yy] = meshgrid(patch_width/2:patch_width:sqrt(Npatch)*patch_width,patch_width/2:patch_width:sqrt(Npatch)*patch_width);
pts = [reshape(xx,Npatch,1),reshape(yy,Npatch,1)];
dists = squareform(pdist(pts));

max_rxns = 1000000;
stoich_matrix = [-1 1;0 -1]; % infection removes an S and adds an I; recovery removes an I
num_species = size(stoich_matrix, 2);

t = cell(num_iter,1);
x = cell(num_iter,1);
XX = cell(num_iter,1);
rxnss = cell(num_iter,1);
rxn_counts = zeros(num_iter,1);

for iter = 1:num_iter
    rxn_count = 1;
    
    T = zeros(max_rxns, 1);
    rxns = zeros(max_rxns,1);
    
    X = zeros(max_rxns,num_species,Npatch);
    X(1,1,:) = s0;
    X(1,2,:) = i0;
    
    f = eye(Npatch); % initial coupling matrix has only intra-patch interactions
    cp = gamma * rep .* n0 ; % contact rate * probability of infection
    
    T(1,:) = tspan(1);
    
    steptau = 5; % time interval between steps
    
    % von mises distribution
    vonMises = @(alpha, thetahat, kappa) exp( kappa*(cos(alpha-thetahat)-1)) / (2*pi*besseli(0,kappa,1));
    
    % expected value of cos theta
    c = @(thetalower, thetaupper) integral(@(theta) cos(theta).*(1/(thetaupper-thetalower)), thetalower, thetaupper);
    
    % expected squared displacement for n steps and constant step size l
    msd = @(n,l,c) n.*l.^2 + 2*l.^2.*(c./(1 - c)).*(n - (1 - c.^n)./(1 - c));
    
    
    %% MAIN LOOP
    while T(rxn_count) < tspan(2) %|| sum(sum(X(rxn_count,:,:))) == 0
        %% Calculate reaction propensities
        beta = repmat(cp,1,Npatch) .* (f ./ repmat(n0',Npatch,1)); % rate of S + I -> 2I
        a = zeros(Npatch, num_species);
        for k = 1:Npatch
            a(k,1) = beta(k,:) * (X(rxn_count,1,k).*squeeze(X(rxn_count,2,:)));
            a(k,2) = gamma*X(rxn_count,2,k);
        end
        a = reshape(a',numel(a),1);
        
        %% Generate two random numbers and calculate timestep tau and reaction mu
        a0 = sum(a);
        tau = -log(rand)/a0; %that is, (1/a0)*log(1/r(1));
        mu = find(cumsum(a) >= rand*a0,1);
        
        if rxn_count + 1 > max_rxns
            t = T(1:rxn_count);
            x = X(1:rxn_count,:,:);
            rxns = rxns(1:rxn_count);
            warning('SSA:ExceededCapacity',...
                'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
            return;
        end
        
        % Update time and carry out reaction mu
        T(rxn_count+1)   = T(rxn_count)   + tau;
        X(rxn_count+1,:,:) = X(rxn_count,:,:); % nothing happened for most of the patches
        X(rxn_count+1,:,ceil(mu/2)) = X(rxn_count,:,ceil(mu/2)) + stoich_matrix(mod(mu,-2)+2,:); % reaction that occurred
        rxns(rxn_count+1) = mu;
        
        %         if mod(rxn_count,5) == 0 % each step is taken at every 10 reactions, though the number of steps is approximated as # reactions * step time interval / actual time
        %             % recalculate coupling matrix
        %             distkernel = @(x) vonMises(x * pi/(max(max(dists))),0,10*2/sqrt(pi*msd(T(rxn_count)/steptau,l,c(-1*thetastar,thetastar)))) * pi/(max(max(dists))) ;
        %             f = distkernel(dists);
        %             % make f row-stochastic?
        %             f = f./repmat(sum(f,2),1,Npatch);
        %         end
        % recalculate coupling matrix
        
        distkernel = @(x) vonMises(x * pi/(max(max(dists))),0,100*2/sqrt(pi*msd(T(rxn_count+1),l,c(-1*thetastar,thetastar)))) * pi/(max(max(dists))) ;
        f = distkernel(dists);
        % make f row-stochastic?
        f = f./repmat(sum(f,2),1,Npatch);
        
        rxn_count = rxn_count + 1;
        
    end
    
    t{iter} = T(1:rxn_count);
    x{iter} = sum(X(1:rxn_count,:,:),3);
    XX{iter} = X(1:rxn_count,:,:);
    rxnss{iter} = rxns(1:rxn_count);
    rxn_counts(iter) = rxn_count;
end


