function [ t, x, y, xx, yy, rxnss, rxnss_simple, rxn_counts, XX, YY ] = gillespie_food_infection_division_iterate( Npatch, s0, i0, sus0, inf0, n0, cp, int_w, p_inf, gamma, l, thetastar, num_iter, patch_width, inf_penalty)
%% Implementation of the Gillespie algorithm for N coupled population SIR-type model for trophallaxis with infection and division of labor
%  using dynamically updated coupling matrix determined by von Mises
%  probability distribution with variance related to mean square
%  displacement
%
%  Populations:
%  Forager (F), Indoor Worker (W), Indoor Young (Y)
%  FS, FI_25, FI_50, FI_75, FI_100, WS, WI_25, WI_50, WI_75, WI_100, YS,
%  YI_25, YI_50, YI_75, YI_100
%  S + I_50 -> I_25 + I_25
%  S + I_75 -> I_25 + I_50
%  S + I_100 -> I_25 + I_75
%  I_25 + I_50 -> I_50 + I_25
%  I_25 + I_75 -> I_50 + I_50
%  I_25 + I_100 -> I_50 + I_75
%  I_50 + I_75 -> I_75 + I_50
%  I_50 + I_100 -> I_75 + I_75
%  I_75 + I_100 -> I_100 + I_75
%  Inf -> R
%
%
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
%       s0:             Initial number hungry (num_groups x Npatch matrix).
%
%       i0:             Initial number full (I_100) (num_groups x Npatch matrix).
%
%       sus0:           Initial number susceptible (num_groups x Npatch matrix).
%
%       inf0:           Initial number infected (num_groups x Npatch matrix).
%
%       n0:             Total population of each patch (num_groups x Npatch matrix).
%
%       cp:             Contact rate * probability of trophallaxis (num_groups x Npatch matrix).
%
%       int_w:          Interaction weights (symmetric num_groups x num_groups matrix).
%
%       p_inf:          Probability of infection occurring from a trophallaxis interaction (num_groups x 1 vector).
%
%       gamma:          Recovery rate from infection
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
num_groups = 3; % F, W, Y
tspan = [0, 10000000];
if nargin < 14
    patch_width = 30; % width of each patch
end
if nargin < 15
    inf_penalty = 1;
end
[xcoord,ycoord] = meshgrid(patch_width/2:patch_width:sqrt(Npatch)*patch_width,patch_width/2:patch_width:sqrt(Npatch)*patch_width);
pts = [reshape(xcoord,Npatch,1),reshape(ycoord,Npatch,1)];
dists = squareform(pdist(pts));

max_rxns = 1000000;
stoich_matrix = [-1  2 -1  0  0;
    -1  1  1 -1  0;
    -1  1  0  1 -1;
    0  0  0  0  0;
    0 -1  2 -1  0;
    0 -1  1  1 -1;
    0  0  0  0  0;
    0  0 -1  2 -1;
    0  0  0  0  0];

stoich_1 = [-1  1  0  0  0;
    -1  1  0  0  0;
    -1  1  0  0  0;
    0 -1  1  0  0;
    0 -1  1  0  0;
    0 -1  1  0  0;
    0  0 -1  1  0;
    0  0 -1  1  0;
    0  0  0 -1  1];

stoich_2 = [ 0  1 -1  0  0;
    0  0  1 -1  0;
    0  0  0  1 -1;
    0  1 -1  0  0;
    0  0  1 -1  0;
    0  0  0  1 -1;
    0  0  1 -1  0;
    0  0  0  1 -1;
    0  0  0  1 -1];

num_species = size(stoich_matrix, 2);
num_rxns = size(stoich_matrix, 1);

t = cell(num_iter,1);
x = cell(num_iter,1);
y = cell(num_iter,1);
xx = cell(num_iter,1);
yy = cell(num_iter,1);
XX = cell(num_iter,1);
YY = cell(num_iter,1);
rxnss = cell(num_iter,1);
rxnss_simple = cell(num_iter,1);
rxn_counts = zeros(num_iter,1);

for iter = 1:num_iter
    iter
    rxn_count = 1;
    
    T = zeros(max_rxns, 1);
    rxns = zeros(max_rxns,1);
    rxns_simple = zeros(max_rxns,1);
    
    X = zeros(max_rxns,num_species,num_groups,Npatch); % tracks fed bees
    X(1,1,:,:) = s0;
    X(1,5,:,:) = i0;
    
    Y = zeros(max_rxns,2,num_groups,Npatch); % tracks infected bees
    Y(1,1,:,:) = sus0;
    Y(1,2,:,:) = inf0;
    
    f = eye(Npatch); % initial coupling matrix has only intra-patch interactions
    %     cp = gamma * rep .* n0 ; % contact rate * probability of infection
    
    T(1,:) = tspan(1);
    
    
    % von mises distribution
    vonMises = @(alpha, thetahat, kappa) exp( kappa*(cos(alpha-thetahat)-1)) / (2*pi*besseli(0,kappa,1));
    
    % expected value of cos theta
    c = @(thetalower, thetaupper) integral(@(theta) cos(theta).*(1/(thetaupper-thetalower)), thetalower, thetaupper);
    
    % expected squared displacement for n steps and constant step size l
    msd = @(n,l,c) n.*l.^2 + 2*l.^2.*(c./(1 - c)).*(n - (1 - c.^n)./(1 - c));
    
    
    %% MAIN LOOP
    while T(rxn_count) < tspan(2) %|| sum(sum(X(rxn_count,:,:))) == 0
        %% Calculate reaction propensities
        beta = permute(repmat(cp,1,1,Npatch),[2 3 1]) .* (f ./ permute(repmat(n0,1,1,Npatch),[2 3 1])); % rate of S + I -> 2I
        beta(isnan(beta)) = 0;
        beta(beta==Inf) = 0;
        a = zeros(Npatch, Npatch, num_groups, num_groups, num_rxns);
        for k = 1:Npatch
            for j = 1:Npatch
                for i = 1:num_groups
                    for h = 1:num_groups
                        a(k,j,i,h,1) = beta(k,j,i) .* X(rxn_count,1,i,k) .* (int_w(i,h).*(X(rxn_count,3,h,j)));
                        a(k,j,i,h,2) = beta(k,j,i) .* X(rxn_count,1,i,k) .* (int_w(i,h).*(X(rxn_count,4,h,j)));
                        a(k,j,i,h,3) = beta(k,j,i) .* X(rxn_count,1,i,k) .* (int_w(i,h).*(X(rxn_count,5,h,j)));
                        a(k,j,i,h,4) = beta(k,j,i) .* X(rxn_count,2,i,k) .* (int_w(i,h).*(X(rxn_count,3,h,j)));
                        a(k,j,i,h,5) = beta(k,j,i) .* X(rxn_count,2,i,k) .* (int_w(i,h).*(X(rxn_count,4,h,j)));
                        a(k,j,i,h,6) = beta(k,j,i) .* X(rxn_count,2,i,k) .* (int_w(i,h).*(X(rxn_count,5,h,j)));
                        a(k,j,i,h,7) = beta(k,j,i) .* X(rxn_count,3,i,k) .* (int_w(i,h).*(X(rxn_count,4,h,j)));
                        a(k,j,i,h,8) = beta(k,j,i) .* X(rxn_count,3,i,k) .* (int_w(i,h).*(X(rxn_count,5,h,j)));
                        a(k,j,i,h,9) = beta(k,j,i) .* X(rxn_count,4,i,k) .* (int_w(i,h).*(X(rxn_count,5,h,j)));
                        if h == 1
                            a(k,j,i,h,:) = a(k,j,i,h,:) * inf_penalty;
                        end
                    end
                end
            end
        end
        a = reshape(permute(a,[2 1 3 4 5]),numel(a),1);
        gammas = zeros(Npatch,num_groups);
        for k = 1:Npatch
            for j = 1:num_groups
                gammas(k,j) = gamma * Y(rxn_count,2,j,k);
            end
        end
        gammas = reshape(gammas,numel(gammas),1);
        a = vertcat(a,gammas);
        
        %% Generate two random numbers and calculate timestep tau and reaction mu
        a0 = sum(a);
        tau = -log(rand)/a0; % that is, (1/a0)*log(1/r(1));
        mu = find(cumsum(a) >= rand*a0,1);
        rxn_no = ceil(mu/(Npatch*Npatch*num_groups*num_groups));
        if rxn_no > num_rxns
            mu_g = mu - Npatch*Npatch*num_groups*num_groups*num_rxns;
            [patchA, groupA] = ind2sub([Npatch,num_groups],mu_g);
        else
            [patchB,patchA,groupA,groupB,~] = ind2sub([Npatch, Npatch, num_groups, num_groups, num_rxns],mu);
        end
        
        %         % this shouldn't really be necessary
        %         if a0 <= 0.0001
        %             break;
        %         end
        %
        %         % this also should not be necessary
        %         if speciesA == 0 | speciesB == 0
        %             break;
        %         end
        
        if rxn_count + 1 > max_rxns
            t = T(1:rxn_count);
            x = X(1:rxn_count,:,:);
            rxns = rxns(1:rxn_count);
            rxns_simple = rxns_simple(1:rxn_count);
            warning('SSA:ExceededCapacity',...
                'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
            return;
        end
        
        % Update time and carry out reaction mu
        T(rxn_count+1)   = T(rxn_count)   + tau;
        X(rxn_count+1,:,:,:) = X(rxn_count,:,:,:); % nothing happened for most of the patches
        Y(rxn_count+1,:,:,:) = Y(rxn_count,:,:,:);
        if rxn_no <= num_rxns
            X(rxn_count+1,:,groupA,patchA) = X(rxn_count,:,groupA,patchA) + stoich_1(rxn_no,:); % reaction that occurred
            X(rxn_count+1,:,groupB,patchB) = X(rxn_count,:,groupB,patchB) + stoich_2(rxn_no,:);
            
            % check if infection occurs
            b = [p_inf(groupA) * Y(rxn_count,1,groupA,patchA) * Y(rxn_count,2,groupB,patchB);
                p_inf(groupB) * Y(rxn_count,1,groupB,patchB) * Y(rxn_count,2,groupA,patchA);
                (1-p_inf(groupA))* Y(rxn_count,1,groupA,patchA) * Y(rxn_count,2,groupB,patchB);
                (1-p_inf(groupB))* Y(rxn_count,1,groupB,patchB) * Y(rxn_count,2,groupA,patchA)];
            mub = find(cumsum(b) >= rand*sum(b),1);
            if b(1) >= 0.00001 | b(2) >= 0.00001
                if mub == 1 % B infects A
                    Y(rxn_count+1,:,groupA,patchA) = Y(rxn_count,:,groupA,patchA) + [-1 1];
                elseif mub == 2 % A infects B
                    Y(rxn_count+1,:,groupB,patchB) = Y(rxn_count,:,groupB,patchB) + [-1 1];
                end
            end
        else
            Y(rxn_count+1,:,groupA,patchA) = Y(rxn_count,:,groupA,patchA) + [0, -1];
        end
        rxns(rxn_count+1) = mu;
        rxns_simple(rxn_count+1) = rxn_no;
        
        distkernel = @(x) vonMises(x * pi/(max(max(dists))),0,100*2/sqrt(pi*msd(T(rxn_count+1),l,c(-1*thetastar,thetastar)))) * pi/(max(max(dists))) ;
        f = distkernel(dists);
        % make f row-stochastic?
        f = f./repmat(sum(f,2),1,Npatch);
        
        rxn_count = rxn_count + 1;
        
    end
    
    t{iter} = T(1:rxn_count);
    x{iter} = sum(X(1:rxn_count,:,:,:),4);
    y{iter} = sum(Y(1:rxn_count,:,:,:),4);
    xx{iter} = sum(sum(X(1:rxn_count,:,:,:),4),3);
    yy{iter} = sum(sum(Y(1:rxn_count,:,:,:),4),3);
    XX{iter} = X(1:rxn_count,:,:,:);
    YY{iter} = Y(1:rxn_count,:,:,:);
    rxnss{iter} = rxns(1:rxn_count);
    rxnss_simple{iter} = rxns_simple(1:rxn_count);
    rxn_counts(iter) = rxn_count;
end




