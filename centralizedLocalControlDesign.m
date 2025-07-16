% This function designs a centralized local control strategy for a microgrid
% comprising multiple Distributed Generators (DGs) and transmission lines.
% The function computes the passivity parameters (P, K, nu, rho, gamma) for
% each DG and line by solving a series of Linear Matrix Inequalities (LMIs).
% This version includes sector-bounded approach for handling CPL nonlinearities.

function [DG,Line,statusLocalController] = centralizedLocalControlDesign(DG,Line,B_il,BarGamma,piVals,plVals)

numOfDGs = size(B_il,1);
numOfLines = size(B_il,2);

epsilon = 0.00001; % Minimum value
useNewCon8 = 0;
debugMode = 0;

% Create LMI variables necessary for (66)
%% Variables corresponding to DGs like P_i, K_i, nu_i, rhoTilde_i, gammaTilde_i
for i = 1:1:numOfDGs
    P_i{i} = sdpvar(3, 3, 'symmetric');
    Kp_i{i} = sdpvar(1, 1, 'full');  % Proportional gain
    Ki_i{i} = sdpvar(1, 1, 'full');  % Integral gain
    nu_i{i} = sdpvar(1, 1, 'full');
    rho_i{i} = sdpvar(1, 1, 'full'); % Representing: rho (for newCon8)
    rhoTilde_i{i} = sdpvar(1, 1, 'full'); % Representing: 1/rho
    gammaTilde_i{i} = sdpvar(1, 1,'full');
    
    % New: For sector-bounded CPL constraints
    lambda_i{i} = sdpvar(1, 1, 'full'); % Lambda parameter for S-procedure
    R_i{i} = sdpvar(3, 3, 'symmetric'); % R matrix for quadratic bound
end

% Variables corresponding to Lines like P_l, nu_l, rho_l
for l = 1:1:numOfLines
    P_l{l} = sdpvar(1, 1, 'symmetric');
    nu_l{l} = sdpvar(1, 1, 'full');
    rho_l{l} = sdpvar(1, 1, 'full'); 
end

constraints = [];

constraintTags = {}; % Cell array to hold tags
constraintMats = {}; % Cell array to hold matrices

%% Combine constraints over all DGs (66a-Part1, 66b,66d,66e)

for i = 1:1:numOfDGs
    K_i{i} = [Kp_i{i}, 0, Ki_i{i}];  % This enforces [Kp, 0, Ki] structure
        
    Ai = DG{i}.A;
    Bi = DG{i}.B;
    Ii = eye(size(Ai));

    % Extract CPL parameters for sector-bounded constraints
    PL_i = DG{i}.PL;     % Constant Power Load value
    Vmin = 0.95 * DG{i}.refVoltage;  % Minimum operating voltage (5% below nominal)
    Vmax = 1.05 * DG{i}.refVoltage;  % Maximum operating voltage (5% above nominal)
    
    % Calculate sector bounds based on Lemma 4
    alpha_i = PL_i / (DG{i}.C * Vmax^2);  % Lower bound
    beta_i = PL_i / (DG{i}.C * Vmin^2);   % Upper bound
    
    % Define the Ti matrix that extracts voltage component (vector form)
    Ti = [1; 0; 0];  % Column vector to extract first element (voltage)
    
    % Create the sector-bounded constraint matrix from Lemma 5
    sector_matrix = [
        -alpha_i*beta_i, (alpha_i+beta_i)/2;
        (alpha_i+beta_i)/2, -1
    ];
    
    % Create the constraint matrices for S-procedure
    % We need to fix the dimension issue in the Kronecker product
    S_matrix = [
        R_i{i}, -P_i{i}*Ti;
        -Ti'*P_i{i}, 0
    ];
    
    % Create the second part of the S-procedure constraint
    S_proc_term = lambda_i{i} * [
        Ti*Ti'*(-alpha_i*beta_i), Ti*((alpha_i+beta_i)/2);
        ((alpha_i+beta_i)/2)*Ti', -1
    ];
    
    % Combine to create the final constraint
    S_constraint = S_matrix - S_proc_term;
    
    slack_nu = sdpvar(1, 1, 'full');
    slack_rho = sdpvar(1, 1, 'full');

    
    tagName = ['gammaTilde_',num2str(i),'_low'];
    constraintTags{end+1} = tagName;
    con0_1 = tag(gammaTilde_i{i} >= epsilon, tagName);
    constraintMats{end+1} = gammaTilde_i{i};

    tagName = ['gammaTilde_',num2str(i),'_high'];
    constraintTags{end+1} = tagName;
    con0_2 = tag(gammaTilde_i{i} <= BarGamma, tagName);
    constraintMats{end+1} = gammaTilde_i{i};

    % Constraint (66a-Part1)
    tagName = ['P_',num2str(i)];
    constraintTags{end+1} = tagName;
    con1 = tag(P_i{i} >= epsilon*Ii, tagName);
    constraintMats{end+1} = P_i{i};
    
    % Constraint for sector-bounded CPL (from Lemma 5)
    tagName = ['S_matrix_',num2str(i)];
    constraintTags{end+1} = tagName;
    con_sector = tag(S_constraint >= epsilon*eye(size(S_constraint)), tagName);
    constraintMats{end+1} = S_constraint;
    
    % Constraint to ensure lambda is positive (for S-procedure)
    tagName = ['lambda_',num2str(i),'_positive'];
    constraintTags{end+1} = tagName;
    con_lambda = tag(lambda_i{i} >= epsilon, tagName);
    constraintMats{end+1} = lambda_i{i};
    
    % Constraint (66b) modified to include R matrix from sector-bounded approach
    DMat = rhoTilde_i{i} * eye(3);
    MMat = [P_i{i}, zeros(3)];
    ThetaMat = [-Ai*P_i{i} - P_i{i}'*Ai' - Bi*K_i{i} - K_i{i}'*Bi' - R_i{i}, -eye(3) + 0.5*P_i{i};
                -eye(3) + 0.5*P_i{i}, -nu_i{i}*eye(3)];
    W = [DMat, MMat; MMat', ThetaMat];
    
    tagName = ['W_',num2str(i)];
    constraintTags{end+1} = tagName;
    con2 = tag(W >= epsilon*eye(size(W)), tagName);
    constraintMats{end+1} = W;
    
    % Constraint (66d)
    p_i{i} = piVals(i); % predefined value

    tagName = ['nu_',num2str(i),'_low'];
    constraintTags{end+1} = tagName;
    con3_1 = tag(nu_i{i} + slack_nu >= -1.5*gammaTilde_i{i}/p_i{i}, tagName);
    constraintMats{end+1} = -gammaTilde_i{i}/p_i{i};

    tagName = ['nu_',num2str(i),'_high'];
    constraintTags{end+1} = tagName;
    con3_2 = tag(nu_i{i} - slack_nu  <= -epsilon/10, tagName);
    constraintMats{end+1} = nu_i{i};

    % Add constraints on slack variables
    tagName = ['slack_nu_',num2str(i)];
    constraintTags{end+1} = tagName;
    con3_1s = tag(slack_nu >= 0, tagName);
    
    costGamma = 0;
    % Add penalty to cost function
    costGamma = costGamma + 1000*slack_nu;
    
    % Constraint (66e)
    tagName = ['rhoTilde_',num2str(i),'_low'];
    constraintTags{end+1} = tagName;
    con4_1 = tag(rhoTilde_i{i} >= epsilon/10, tagName);
    constraintMats{end+1} = rhoTilde_i{i};

    tagName = ['rhoTilde_',num2str(i),'_high1'];
    constraintTags{end+1} = tagName;
    con4_21 = tag(rhoTilde_i{i} <= 1.5*p_i{i}, tagName);
    constraintMats{end+1} = p_i{i};
    
    tagName = ['rhoTilde_',num2str(i),'_high2'];
    constraintTags{end+1} = tagName;
    con4_22 = tag(rhoTilde_i{i} <= 4*gammaTilde_i{i}/p_i{i}, tagName);
    constraintMats{end+1} = 4*gammaTilde_i{i}/p_i{i};

    % New con8:
    if useNewCon8
        tagName = ['rho_',num2str(i),'_low'];
        constraintTags{end+1} = tagName;
        con4_31 = tag(rho_i{i} >= 1/max(p_i{i}, 4*BarGamma/p_i{i}), tagName);
        constraintMats{end+1} = 1/max(p_i{i}, 4*BarGamma/p_i{i});
    
        tagName = ['rho_',num2str(i),'_high'];
        constraintTags{end+1} = tagName;
        con4_32 = tag(rho_i{i} <= 1/epsilon, tagName);
        constraintMats{end+1} = rho_i{i};
    
        tagName = ['rho_',num2str(i),'_receprocity'];
        constraintTags{end+1} = tagName;
        con4_33 = tag([rho_i{i}, 1; 1, rhoTilde_i{i}] >= 0, tagName); % rho_i rhoTilde_i >= 1
        constraintMats{end+1} = [rho_i{i}, 1; 1, rhoTilde_i{i}];
    end

    % Collecting Constraints
    if useNewCon8
        constraints = [constraints, con0_1, con0_2, con1, con2, con3_1, con3_2, con3_1s, con4_1, con4_21, con4_22, con4_31, con4_32, con4_33, con_sector, con_lambda];
    else
        constraints = [constraints, con0_1, con0_2, con1, con2, con3_1, con3_2, con3_1s, con4_1, con4_21, con4_22, con_sector, con_lambda];
    end
end

%% Combine constraints over all Lines (66a-Part2, 66c)
for l = 1:1:numOfLines
    Rl = Line{l}.R;
    Ll = Line{l}.L;
    Il = eye(1);
    % Constraint (66a-Part2)

    tagName = ['PBar_',num2str(l)];
    constraintTags{end+1} = tagName;
    con5 = tag(P_l{l} >= epsilon*Il, tagName);
    constraintMats{end+1} = P_l{l};
    
    p_l{l} = plVals(l); %1/numOfLines;  % predefined value

    
    % Constraint (66c)
    W = [(2*P_l{l}*Rl)/Ll - rho_l{l}, -P_l{l}/Ll + 1/2;
         -P_l{l}/Ll + 1/2, -nu_l{l}];

    tagName = ['WBar_',num2str(l)];
    constraintTags{end+1} = tagName;
    con6 = tag(W >= epsilon*eye(size(W)), tagName);
    constraintMats{end+1} = W;
    
    tagName = ['nuBar_',num2str(l),'_high'];
    constraintTags{end+1} = tagName;
    con7_0 = tag(nu_l{l} <= -epsilon, tagName);
    constraintMats{end+1} = nu_l{l};

    tagName = ['rhoBar_',num2str(l),'low'];
    constraintTags{end+1} = tagName;
    con7_1 = tag(rho_l{l} >= epsilon, tagName);
    constraintMats{end+1} = rho_l{l};

    % Collecting Constraints
    constraints = [constraints, con5, con6, con7_0, con7_1];
end

%% Combine all mixed constraints (66f, 66g)
for i = 1:1:numOfDGs
    for l = 1:1:numOfLines
        Ct = DG{i}.C;

        if B_il(i,l) ~= 0
            slack1 = sdpvar(1, 1, 'full');
            slack2 = sdpvar(1, 1, 'full');
       
           % Constraint (66f)
           tagName = ['rhoBar_',num2str(l),'_low1_',num2str(i)];
           constraintTags{end+1} = tagName;
           con7_2 = tag(rho_l{l} + slack1 >= -(p_i{i}*nu_i{i})/(p_l{l}*Ct^2), tagName);
           constraintMats{end+1} = -(p_i{i}*nu_i{i})/(p_l{l}*Ct^2);
           
           tagName = ['rhoBar_',num2str(l),'_low2_',num2str(i)];
           constraintTags{end+1} = tagName;
           con7_3 = tag(rho_l{l} + slack2 >= ((rhoTilde_i{i})/(p_i{i}*p_l{l}))*((p_i{i}/(2*Ct))-(p_l{l}/2))^2, tagName);
           constraintMats{end+1} = ((rhoTilde_i{i})/(p_i{i}*p_l{l}))*((p_i{i}/(2*Ct))-(p_l{l}/2))^2;

           % Add constraints on slack variables
           tagName = ['slack1_',num2str(l),'_',num2str(i)];
           constraintTags{end+1} = tagName;
           con7_2s = tag(slack1 >= 0, tagName);
           
           tagName = ['slack2_',num2str(l),'_',num2str(i)];
           constraintTags{end+1} = tagName;
           con7_3s = tag(slack2 >= 0, tagName);
                
           constraints = [constraints, con7_2, con7_3, con7_3s, con7_2s];

           costGamma = costGamma + 1000*(slack1 + slack2);

           if useNewCon8
               % New con8 (66g):
               tagName = ['nuBar_',num2str(l),'_low_',num2str(i)];
               constraintTags{end+1} = tagName;
               con8New = nu_l{l} >= -p_i{i}*rho_i{i}/p_l{l};
               constraintMats{end+1} = -p_i{i}*rho_i{i}/p_l{l};
           else
               % Constraint (66g)
               n = 1;          % Number of intervals
               rho_max = min(p_i{i}, 4*BarGamma/p_i{i});
               rho_min = rho_max/1000;
    
               delta_i = (rho_max - rho_min) / n;
               
               % Initialize cell array to store individual constraints
               con8 = [];
    
               % Loop over each k from 1 to n to create individual constraints
               tilde_rho_i_prev = rho_min;
               tilde_y_i_prev = -p_i{i} / (p_l{l} * tilde_rho_i_prev);
                
               for k = 1:n
                   % Compute tilde_rho_i^k
                   tilde_rho_i_k = rho_min + k * delta_i;
    
                   % Compute tilde_y_i^k
                   tilde_y_i_k = -p_i{i} / (p_l{l} * tilde_rho_i_k);
    
                   % Compute m_k and c_k
                   m_k = (tilde_y_i_k - tilde_y_i_prev) / delta_i;
                   c_k = tilde_y_i_k - m_k * tilde_rho_i_k;
    
                   % Define Constraint (66g)
                   tagName = ['nuBar_',num2str(l),'_low',num2str(n),'_',num2str(i)];
                   constraintTags{end+1} = tagName;
                   con8_k = tag(nu_l{l} >= m_k * rhoTilde_i{i} + c_k, tagName);
                   constraintMats{end+1} = m_k * rhoTilde_i{i} + c_k;
    
                   con8 = [con8, con8_k];
                   
                   % Compute tilde_rho_i^{k-1} and tilde_y_i^{k-1}
                   tilde_rho_i_prev = tilde_rho_i_k;
                   tilde_y_i_prev = tilde_y_i_k;
               end
           end

           if useNewCon8
               constraints = [constraints, con8New];
           else
               constraints = [constraints, con8];
           end
        end
    end
end

%% Cost function formulation
costGamma = 0;
for i = 1:numOfDGs
    if useNewCon8
        costGamma = costGamma + (-10*nu_i{i}+rhoTilde_i{i}) + gammaTilde_i{i} + trace(P_i{i}) - 1000000*rho_i{i};
    else
        costGamma = costGamma + (-10*nu_i{i}+rhoTilde_i{i}) + gammaTilde_i{i} - 100*trace(P_i{i});
    end
    
    % Add cost term for sector-bounded constraints
    costGamma = costGamma + 100*lambda_i{i} + 10*trace(R_i{i});
end

for l = 1:numOfLines
    costGamma = costGamma + (-100*nu_l{l}-rho_l{l}) + 1*trace(P_l{l});
end

% Defining costfunction
costFunction = costGamma;

% Add a regularization term to penalize large controller gains
for i = 1:numOfDGs
    costFunction = costFunction + 1 * norm(K_i{i}, 'fro');  % Frobenius norm of the controller gain
end

solverOptions = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug', 0);

sol = optimize(constraints, costFunction, solverOptions);

statusLocalController = sol.problem == 0;

%% Extract variable values
for i = 1:1:numOfDGs
    P_iVal = value(P_i{i});
    Kp_iVal = value(Kp_i{i});
    Ki_iVal = value(Ki_i{i});
    nu_iVal = value(nu_i{i});
    rhoTilde_iVal = value(rhoTilde_i{i});
    rho_iVal = 1 / rhoTilde_iVal;
    gammaTilde_iVal = value(gammaTilde_i{i});
    R_iVal = value(R_i{i});  % Get R matrix for CPL constraint
    lambda_iVal = value(lambda_i{i});  % Get lambda for S-procedure

    Kp0 = Kp_iVal / P_iVal(1,1);  % Use only the (1,1) element of P for Kp
    Ki0 = Ki_iVal / P_iVal(3,3);  % Use only the (3,3) element of P for Ki
        
    % Construct K0 with guaranteed zero middle element
    K0 = [Kp0, 0, Ki0];

    % update DG
    DG{i}.P0 = P_iVal;
    DG{i}.K0 = K0;
    DG{i}.nu = nu_iVal;
    DG{i}.rho = rho_iVal;
    DG{i}.gammaTilde0 = gammaTilde_iVal;
    DG{i}.R_CPL = R_iVal;  % Store R matrix for CPL constraint
    DG{i}.lambda_CPL = lambda_iVal;  % Store lambda for S-procedure
    
    % Also store the sector bounds for CPL
    PL_i = DG{i}.PL;
    Vmin = 0.95 * DG{i}.refVoltage;
    Vmax = 1.05 * DG{i}.refVoltage;
    DG{i}.alpha_CPL = PL_i / (DG{i}.C * Vmax^2);
    DG{i}.beta_CPL = PL_i / (DG{i}.C * Vmin^2);
end

for l = 1:1:numOfLines
    P_lVal = value(P_l{l});
    nu_lVal = value(nu_l{l});
    rho_lVal = value(rho_l{l});

    % update Line
    Line{l}.P0 = P_lVal;
    Line{l}.nu = nu_lVal;
    Line{l}.rho = rho_lVal;
end

% Debug mode for checking violated constraints
if debugMode
    feasibility = check(constraints);
    
    % Combine tags and feasibility into one array for sorting
    combinedList = [feasibility(:), (1:length(feasibility))'];
    
    % Sort based on the first column (feasibility values)
    sortedList = sortrows(combinedList, 1);  % Sort by feasibility, ascending order
    
    % Printing
    for i = 1:length(sortedList)
        idx = sortedList(i, 2);  % Get the original index of the constraint
        if feasibility(idx) < -1e-6
            disp(['Constraint "', constraintTags{idx}, '" is violated by ', num2str(feasibility(idx)), ' .']);
            W_val = value(constraintMats{idx});
            for j = 1:size(W_val, 1)
                submatrix = W_val(1:j, 1:j);  % Extract principal minor
                if det(submatrix) < 0
                    disp(['Principal minor ', num2str(j), ' is not positive semi-definite.']);
                end
            end
        else
            disp(['Constraint "', constraintTags{idx}, '" is satisfied by ',num2str(feasibility(idx)),' .']);
            W_val = value(constraintMats{idx});
        end
    end

    % Checking the nonlinear constraint
    for i = 1:1:numOfDGs
        for l = 1:1:numOfLines
            Ct = DG{i}.C;
            if B_il(i,l) ~= 0
                nonLinCons = value(nu_l{l} +  p_i{i}/(rhoTilde_i{i}*p_l{l}));
                del3 = value((rho_l{l}*Ct^2/(-nu_i{i})) - ((-nu_l{l})*rhoTilde_i{i}));
                
                T0 = Ct*(Ct*rho_l{l}/rhoTilde_i{i} + 2);
                T1 = (-nu_l{l})*rhoTilde_i{i};
                T2 = (-nu_i{i})/rho_l{l};
                del4 = value(T0-(T1+T2));
                
                % Also check CPL sector bounds
                PL_i = DG{i}.PL;
                Vmin = 0.95 * DG{i}.refVoltage;
                Vmax = 1.05 * DG{i}.refVoltage;
                alpha_i = PL_i / (DG{i}.C * Vmax^2);
                beta_i = PL_i / (DG{i}.C * Vmin^2);
                
                disp(['DG ', num2str(i), ' CPL sector bounds: alpha=', num2str(alpha_i), ', beta=', num2str(beta_i)]);
            end
        end
    end
end
end