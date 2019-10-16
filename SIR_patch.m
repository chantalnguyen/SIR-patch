%% Solve SIR model using IE method for arbitrary number of patches
% The simulation stops at the end of the epidemic, which is defined to
% occur when the mean number of infectives reaches 0.
%
% Specify:
% Npatch:       number of patches
% s0:           initial numbers of susceptible individuals (Npatch x 1 vector)
% i0:           initial numbers of infected individuals (Npatch x 1 vector)
% rep:          reproductive number of each patch (Npatch x 1 vector)
% f:        	Npatch x Npatch symmetric matrix specifying fractions of home
%               (diagonal terms) and foreign (off-diagonal) contacts. each
%               row & column must sum to 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%      Exact numerical integration of the master     %%%%%%%%%%%%
%%%%%%%%%%%       equation in some models of stochastic        %%%%%%%%%%%%
%%%%%%%%%%%                     epidemiology                   %%%%%%%%%%%%
%%%%%%%%%%%                      IE METHOD                     %%%%%%%%%%%%
%%%%%%%%%%%        Original authors: Jenkinson/Goutsias        %%%%%%%%%%%%
%%%%%%%%%%%             Modified by: Chantal Nguyen            %%%%%%%%%%%%
%%%%%%%%%%%              Last Modified: 10/07/2019             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

Npatch = 5; % number of patches
s0 = [3; 2; 3; 1; 3]; % number of initial susceptible for each patch
i0 = [0; 1; 0; 2; 0]; % number of initial infected for each patch
rep = 10 * ones(Npatch,1);
f = [0, .0864, .0237, .004, .0001;
    .0864, 0, .0864, .0237, .004;
    .0237, .0864, 0, .0864, .0237;
    .004, .0237, .0864, 0, .0864;
    .0001, .004, .0237, .0864, 0]; % coupling matrix
f = f(1:Npatch,1:Npatch);
f(logical(eye(Npatch))) = 1 - sum(f,2); % set diagonals so that matrix is stochastic

tau = 0.1; % timestep in days
tf = 250; % simulation time in days

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
r0 = zeros(Npatch,1); % initial number of recovered
n0 = s0+i0+r0; % total population (constant)

num_rxns = Npatch * 2; % number of reactions (infection and recovery)

gamma = 0.1; % 1/day, rate of I -> R
cp = gamma * rep .* n0 ; % contact rate * probability of infection
beta = repmat(cp,1,Npatch) .* (f ./ repmat(n0',Npatch,1)); % rate of S + I -> 2I

steps = tf/tau; % number of steps

% determine sample space Z
hyperCubeOrigin = zeros(num_rxns,1);       % origin of hypercube
hyperCubeExtent = [s0;s0+i0];              % extent of hypercube = [ max # infections; max # recoveries]
K               = prod(hyperCubeExtent+1); % number of states

% define initial conditions for probability vector P
pPrev    = zeros(K,1);
pPrev(1) = 1;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% create genMatrix file if it does not exist
if ~isfile(['genMatrix_' num2str(Npatch) '.m'])
    fid = fopen(['genMatrix_' num2str(Npatch) '.m'],'w');
    fprintf(fid,'function IminusTauA = genMatrix');
    fprintf(fid,['_' num2str(Npatch)]);
    fprintf(fid,'(beta,gamma,s0,i0,K,hyperCubeOrigin,hyperCubeExtent,num_rxns,tau,propensity)\n');
    fprintf(fid, '%% find hypercube edge lengths\nedge = 1+hyperCubeExtent-hyperCubeOrigin;\n');
    fprintf(fid, '%% initialize\nnonZeroTotal    = K*(num_rxns+1); \n');
    fprintf(fid, 'nonZeroElements = zeros(nonZeroTotal,1);\nrows            = nonZeroElements;\ncols            = rows;\n\n');
    fprintf(fid, 'nonZeroCount = 1;\ncount = 1;\n\n');
    for i = 1:Npatch*2
        fprintf(fid,['for z' num2str(i) '=hyperCubeOrigin(' num2str(i) '):hyperCubeExtent(' num2str(i) ')\n']);
    end
    fprintf(fid,'%% compute current propensities\nalpha=propensity(');
    for i = 1:Npatch*2
        fprintf(fid,['z' num2str(i) ',']);
    end
    fprintf(fid,'beta,gamma,s0,i0);\n');
    fprintf(fid,['%% find diagonal term\nrows(nonZeroCount)            = count;\ncols(nonZeroCount)            = count;\n'...
        'nonZeroElements(nonZeroCount) = tau*sum(alpha);\nnonZeroCount                  = nonZeroCount+1;\n\n']);
    for i = 1:(Npatch*2-1)
        fprintf(fid,['%% find m=' num2str(i) ' off diagonal\n' ...
            'if alpha(' num2str(i) ')>0\n' ...
            'rows(nonZeroCount)            = count+']);
        for j = (i+1):(Npatch*2-1)
            fprintf(fid,['edge(' num2str(j) ')*']);
        end
        fprintf(fid,['edge(' num2str(Npatch*2) ');\n']);
        fprintf(fid,['cols(nonZeroCount)            = count;\n'...
            'nonZeroElements(nonZeroCount) = -tau*alpha(' num2str(i) ');\n'...
            'nonZeroCount                  = nonZeroCount+1;\nend\n']);
    end
    fprintf(fid,['%% find m=' num2str(Npatch*2) ' off diagonal\n']);
    fprintf(fid,['if alpha(' num2str(Npatch*2) ')>0\n'...
        'rows(nonZeroCount)            = count+1;\n']);
    fprintf(fid,['cols(nonZeroCount)            = count;\n'...
        'nonZeroElements(nonZeroCount) = -tau*alpha(' num2str(Npatch*2) ');\n'...
        'nonZeroCount                  = nonZeroCount+1;\nend\n']);
    
    fprintf(fid,'count=count+1;\n\n');
    for i = 1:Npatch*2
        fprintf(fid,'end\n');
    end
    fprintf(fid,['%% generate sparse matrix\n'...
        'minusTauA  = sparse(rows(nonZeroElements~=0),cols(nonZeroElements~=0),...\n'...
        'nonZeroElements(nonZeroElements~=0),K,K);\n'...
        '%% add to identity matrix\n'...
        'IminusTauA = minusTauA + eye(K,''like'',minusTauA);\n']);
    
    fclose(fid);
    
end

% create propensity file if it does not exist
if ~isfile(['propensity_' num2str(Npatch) '.m'])
    fid = fopen(['propensity_' num2str(Npatch) '.m'],'w');
    fprintf(fid,'function alpha = propensity');
    fprintf(fid,['_' num2str(Npatch) '(']);
    for i = 1:Npatch*2
        fprintf(fid,['z' num2str(i) ',']);
    end
    fprintf(fid, 'beta,gamma,s0,i0)\n\n');
    fprintf(fid, ['alpha = zeros(' num2str(Npatch*2) ',1);\n\n']);
    for i = 1:Npatch
        fprintf(fid,['s' num2str(i) ' = s0(' num2str(i) ') - z' num2str(i) ';\n']);
        fprintf(fid,['i' num2str(i) ' = i0(' num2str(i) ') + z' num2str(i) ' - z' num2str(i+Npatch) ';\n']);
    end
    fprintf(fid, '\nss = [');
    for i = 1:Npatch-1
        fprintf(fid, ['s' num2str(i) ';']);
    end
    fprintf(fid,['s' num2str(Npatch) '];\n']);
    fprintf(fid, 'ii = [');
    for i = 1:Npatch-1
        fprintf(fid, ['i' num2str(i) ';']);
    end
    fprintf(fid,['i' num2str(Npatch) '];\n\n']);
    fprintf(fid, 'if ');
    for i = 1:(Npatch-1)
        fprintf(fid, ['z' num2str(i+Npatch) ' <= i0(' num2str(i) ') + z' num2str(i) ' && ']);
    end
    fprintf(fid, ['z' num2str(Npatch*2) ' <= i0(' num2str(Npatch) ') + z' num2str(Npatch) '\n']);
    
    for i = 1:(Npatch)
        fprintf(fid, ['alpha(' num2str(i) ') = beta(' num2str(i) ',:) * (ss(' num2str(i) ').*ii);\n']);
    end
    fprintf(fid, ['alpha(' num2str(Npatch+1) ':' num2str(Npatch*2) ') = gamma.*ii;\n']);
    fprintf(fid, 'else \n');
    fprintf(fid,['alpha(1:' num2str(Npatch*2) ') = 0;\n']);
    fprintf(fid, 'end');
    fclose(fid);
end

% create mapProbs file if it does not exist
% I'm using a different lexicographic ordering of the probability vector
% than that of the original script, which results in a backwards (or
% forwards, I suppose) ordering of the for loops.
% This code takes the K x 1 probability vector P and maps it to a (s0(i)+
% 1) x (s0(i)+i0(i)+1) matrix for each patch i, compiled into an Npatch x 1
% cell. For each patch, the matrix represents the probabilities of being
% in every possible (S,I) state. Columns represent probabilities of having
% each possible number of I, rows represent each possible S. Numbering
% obeys the MATLAB shift. For instance, the sum of the first column
% represents the probability of having 0 infectives.

if ~isfile(['mapProbs_' num2str(Npatch) '.m'])
    fid = fopen(['mapProbs_' num2str(Npatch) '.m'],'w');
    fprintf(fid,['function SIprobs = mapProbs_' num2str(Npatch) '(probVect,hyperCubeOrigin,hyperCubeExtent,s0,i0)\n\n']);
    for j = 1:Npatch
        fprintf(fid,['SIprobs' num2str(j) ' = zeros(s0(' num2str(j) ')+1,s0(' num2str(j) ')+i0(' num2str(j) ')+1);\n']);
    end
    fprintf(fid,'count=1;\n\n');
    for j = 1:Npatch
        for i = 1:Npatch*2
            fprintf(fid,['for z' num2str(i) '=hyperCubeOrigin(' num2str(i) '):hyperCubeExtent(' num2str(i) ')\n']);
        end
        fprintf(fid, ['if probVect(count) > 0 && 1+s0(' num2str(j) ')-z' num2str(j) '>0 && 1+i0(' num2str(j) ')+z' num2str(j) '-z' num2str(j+Npatch) '> 0\n']);
        fprintf(fid, ['SIprobs' num2str(j) '(1+s0(' num2str(j) ')-z' num2str(j) ',1+i0(' num2str(j) ')+z' num2str(j) '-z' num2str(j+Npatch) ') ...\n']);
        fprintf(fid, ['= SIprobs' num2str(j) '(1+s0(' num2str(j) ')-z' num2str(j) ',1+i0(' num2str(j) ')+z' num2str(j) '-z' num2str(j+Npatch) ') ...\n']);
        fprintf(fid, '+ probVect(count);\nend\ncount=count+1;\n');
        for i = 1:Npatch*2
            fprintf(fid,'end\n');
        end
        fprintf(fid,'count = 1;\n');
    end
    fprintf(fid,'SIprobs = {');
    for j = 1:Npatch-1
        fprintf(fid,['SIprobs' num2str(j) ';']);
    end
    fprintf(fid,['SIprobs' num2str(Npatch) '};\n']);
    fclose(fid);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tic;
genMatrix = str2func(['genMatrix_' num2str(Npatch)]);
propensity = str2func(['propensity_' num2str(Npatch)]);
mapProbs = str2func(['mapProbs_' num2str(Npatch)]);

meanS = zeros(steps, Npatch);
meanI = zeros(steps, Npatch);
stdevS = zeros(steps, Npatch);

% find implicit Euler matrix
IminusTauA = genMatrix(beta,gamma,s0,i0,K,hyperCubeOrigin,hyperCubeExtent,num_rxns,tau,propensity);

% perform implicit Euler solver steps
for index = 1:steps
    
    pNext = IminusTauA\pPrev; % solves IminusTauA*pNext = pPrev for pNext
    pPrev = pNext;            % prepare for next time step
    % calculate mean S and I values
    SIprobs = mapProbs(pNext,hyperCubeOrigin,hyperCubeExtent,s0,i0);
    
    for i = 1:Npatch
        meanS(index,i) = (0:s0(i))*sum(SIprobs{i},2);
        stdevS(index,i) = sqrt((((0:s0(i))-meanS(index,i)).^2)*sum(SIprobs{i},2));
        meanI(index,i) = sum(SIprobs{i},1)*(0:(s0(i)+i0(i)))';
    end
    
    % stop simulation when epidemic has ended
    endcrit = 0;
    for i = 1:Npatch
        endcrit = endcrit + sum(SIprobs{i}(1:end,1));
    end
    if endcrit > 0.9999*Npatch
        endtime = index*tau;
        disp(['epidemic has ended after ' num2str(endtime) ' days']);
        break;
    end
    
end

toc;



