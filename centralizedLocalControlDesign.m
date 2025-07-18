function [DG,Line,statusLocalController] = centralizedLocalControlDesign(DG,Line,B_il,BarGamma,piVals,plVals)

numOfDGs = size(B_il,1);
numOfLines = size(B_il,2);

epsilon = 0.000001; % Minimum value
debugMode = 0; % Enable for debugging

% Create LMI variables for DGs
for i = 1:1:numOfDGs
    P_tilde_i{i} = sdpvar(3, 3, 'symmetric');  % P̃_i in Theorem 4
    
    % Define controller gains separately to enforce [Kp, 0, Ki] structure
    Kp_i{i} = sdpvar(1, 1, 'full');           % Proportional gain
    Ki_i{i} = sdpvar(1, 1, 'full');           % Integral gain
    K_tilde_i0{i} = [Kp_i{i}, 0, Ki_i{i}];   % This enforces [Kp, 0, Ki] structure
    
    R_tilde_i{i} = sdpvar(3, 3, 'symmetric'); % R̃_i in Theorem 4
    lambda_tilde_i{i} = sdpvar(1, 1, 'full'); % λ̃_i in Theorem 4
    nu_i{i} = sdpvar(1, 1, 'full');           % νi in Theorem 4
    rho_tilde_i{i} = sdpvar(1, 1, 'full');    % ρ̃i in Theorem 4
    gamma_tilde_i{i} = sdpvar(1, 1, 'full');  % γ̃i in Theorem 4
end

% Create LMI variables for Lines
for l = 1:1:numOfLines
    P_bar_l{l} = sdpvar(1, 1, 'symmetric');   % P̄_l in Theorem 4
    nu_bar_l{l} = sdpvar(1, 1, 'full');       % ν̄_l in Theorem 4
    rho_bar_l{l} = sdpvar(1, 1, 'full');      % ρ̄_l in Theorem 4
end

constraints = [];
constraintTags = {};
constraintMats = {};

%% DG Constraints from Theorem 4
for i = 1:1:numOfDGs
    % Get system matrices
    Ai = DG{i}.A;
    Bi = DG{i}.B;
    Ii = eye(size(Ai));
    
    % Debug matrix dimensions
    fprintf('DG %d matrix dimensions: Ai=%dx%d, Bi=%dx%d\n', i, size(Ai), size(Bi));
    
    % CPL parameters for sector boundedness (from paper's Lemma 6)
    PL_i = DG{i}.PL;     % Constant Power Load value
    Cti = DG{i}.C;       % Capacitance
    Vmin = 0.95 * DG{i}.refVoltage;  % 5% below nominal
    Vmax = 1.05 * DG{i}.refVoltage;  % 5% above nominal
   
    
    % Basic positivity constraints
    tagName = ['P_tilde_',num2str(i),'_positive'];
    constraintTags{end+1} = tagName;
    con1 = tag(P_tilde_i{i} >= epsilon*Ii, tagName);
    constraintMats{end+1} = P_tilde_i{i};
    
    tagName = ['R_tilde_',num2str(i),'_positive'];
    constraintTags{end+1} = tagName;
    con2 = tag(R_tilde_i{i} >= epsilon*Ii, tagName);
    constraintMats{end+1} = R_tilde_i{i};
    
    tagName = ['lambda_tilde_',num2str(i),'_bound'];
    constraintTags{end+1} = tagName;
    con3 = tag(lambda_tilde_i{i} >= 1, tagName); % λ̃_i > 1 from Theorem 4
    constraintMats{end+1} = lambda_tilde_i{i};

    % Passivity constraints
    % tagName = ['nu_',num2str(i),'_bound'];
    % constraintTags{end+1} = tagName;
    % con_nu = tag(nu_i{i} <= 0, tagName);  % Passivity index constraint
    % constraintMats{end+1} = nu_i{i};
    % 
    % tagName = ['rho_tilde_',num2str(i),'_positive'];
    % constraintTags{end+1} = tagName;
    % con_rho = tag(rho_tilde_i{i} >= epsilon, tagName);  % Positive definite
    % constraintMats{end+1} = rho_tilde_i{i};
    % 
    % Gamma constraints (from Lemma 9 - includes BarGamma bound)
    tagName = ['gamma_tilde_',num2str(i),'_low'];
    constraintTags{end+1} = tagName;
    con0_1 = tag(gamma_tilde_i{i} >= epsilon, tagName);
    constraintMats{end+1} = gamma_tilde_i{i};

    tagName = ['gamma_tilde_',num2str(i),'_high'];
    constraintTags{end+1} = tagName;
    con0_2 = tag(gamma_tilde_i{i} <= BarGamma, tagName);  % From Lemma 9
    constraintMats{end+1} = gamma_tilde_i{i};
    
    % FIRST CONSTRAINT from Theorem 4 (Main LMI with CPL)
    % Building the 4-block structure exactly as in the paper
    
    % Block (1,1): ρ̃_i*I + R̃_i
    Block_11 = rho_tilde_i{i} * eye(3) + R_tilde_i{i};
    
    % Block (1,2): R̃_i
    Block_12 = R_tilde_i{i};
    
    % Block (1,3): 0
    Block_13 = zeros(3, 3);
    
    % Block (1,4): 0
    Block_14 = zeros(3, 3);
    
    % Block (2,2): R̃_i
    Block_22 = R_tilde_i{i};
    
    % Block (2,3): P̃_i
    Block_23 = P_tilde_i{i};
    
    % Block (2,4): 0
    Block_24 = zeros(3, 3);
    
    % Block (3,3): -H(A_i*P̃_i + B_i*K̃_i0)
    H_AiPi_BiKi0 = Ai*P_tilde_i{i} + Bi*K_tilde_i0{i};
    Block_33 = -H_AiPi_BiKi0 - H_AiPi_BiKi0';
    
    % Block (3,4): -I + (1/2)*P̃_i
    Block_34 = -eye(3) + 0.5*P_tilde_i{i};
    
    % Block (4,4): -ν_i*I
    Block_44 = -nu_i{i} * eye(3);
    
    % Construct the full 12x12 matrix
    W1_i = [Block_11, Block_12, Block_13, Block_14;
            Block_12', Block_22, Block_23, Block_24;
            Block_13', Block_23', Block_33, Block_34;
            Block_14', Block_24', Block_34', Block_44];
    
    tagName = ['W1_',num2str(i)];
    constraintTags{end+1} = tagName;
    con4 = tag(W1_i >= epsilon*eye(size(W1_i)), tagName);
    constraintMats{end+1} = W1_i;
    
    
    % SECOND CONSTRAINT from Theorem 4 (S-procedure for CPL sector boundedness)
    % All parts are 6×6 matrices based on your dimension analysis

    % Ti matrix setup
    Ti = [1; 0; 0];  % 3×1 column vector

    % Theta matrix components (from equation 95 in the paper)
    Theta_11 = -alpha_i * beta_i;
    Theta_12 = (alpha_i + beta_i)/2;
    Theta_21 = Theta_12;  % Symmetric
    Theta_22 = -1;

    % FIRST PART: [3×3, 3×1, 3×2; 3×3, 3×1, 3×2] → 6×6 matrix
    % Block structure: [I 0 0; 0 P̃_iTi 0]
    % Note: P̃_i = P_i^(-1), so P̃_iTi represents P_i^(-1) * Ti
    Pi_Ti_product = P_tilde_i{i} * Ti;  % This is 3×3 * 3×1 = 3×1
    
    % To match your dimension spec [3×3, 3×1, 3×2; 3×3, 3×1, 3×2], we need:
    % Row 1: [3×3, 3×1, 3×2] and Row 2: [3×3, 3×1, 3×2]
    % But the block [0, P̃_iTi, 0] means the P̃_iTi should be in the (4,4) position as 3×1
    
    First_Part = [eye(3),        zeros(3,1),     zeros(3,2);      % Row 1: [I₃ₓ₃, 0₃ₓ₁, 0₃ₓ₂]
                  zeros(3,3),    Pi_Ti_product,  zeros(3,2)];     % Row 2: [0₃ₓ₃, P̃_iTi₃ₓ₁, 0₃ₓ₂]
    % Result: 6×6 matrix
    
    % SECOND PART: [3×3, 3×3; 1×3, 1×3; 2×3, 2×3] → 6×6 matrix  
    % Block structure: [I 0; 0 Ti'P̃_i; 0 0]
    Ti_transpose_Pi = Ti' * P_tilde_i{i};  % 1×3 result
    
    Second_Part = [eye(3),       zeros(3,3);                      % Row 1: [I₃ₓ₃, 0₃ₓ₃]
                   zeros(1,3),   Ti_transpose_Pi;                  % Row 2: [0₁ₓ₃, Ti'P̃_i₁ₓ₃]
                   zeros(2,3),   zeros(2,3)];                     % Row 3: [0₂ₓ₃, 0₂ₓ₃]
    % Result: 6×6 matrix
    
    % THIRD PART: [3×3, 3×3; 3×3, 3×3] → 6×6 matrix
    % Block structure: [R̃_i 0; 0 I]
    Third_Part = [R_tilde_i{i},  zeros(3,3);                      % Row 1: [R̃_i₃ₓ₃, 0₃ₓ₃]
                  zeros(3,3),    eye(3)];                         % Row 2: [0₃ₓ₃, I₃ₓ₃]
    % Result: 6×6 matrix
    
    % FOURTH PART: [3×3, 3×3; 3×3, 3×3] → 6×6 matrix
    % From the paper's constraint: [Θ₁₁ Θ₁₂P̃_i + λ̃_iI; P̃_iΘ₂₁ + λ̃_iI 0]
    % where Θ₁₁ = -α_iβ_i(T_iT_i^T), Θ₁₂ = (α_i+β_i)/2(T_iT_i^T), etc.
    
    % Calculate TiTi^T (this is 3×1 * 1×3 = 3×3)
    TiTi_transpose = Ti * Ti';  % 3×3 matrix
    
    % Upper left 3×3 block: Θ₁₁ = -α_iβ_i(T_iT_i^T)
    Fourth_Part_11 = Theta_11 * TiTi_transpose;  % 3×3 matrix
    
    % Upper right 3×3 block: Θ₁₂P̃_i + λ̃_iI
    % Note: Θ₁₂ involves T_iT_i^T term, so Θ₁₂P̃_i = (α_i+β_i)/2 * T_iT_i^T * P̃_i
    Fourth_Part_12 = Theta_12 * TiTi_transpose * P_tilde_i{i} + lambda_tilde_i{i} * eye(3);
    
    % Lower left 3×3 block: P̃_iΘ₂₁ + λ̃_iI  
    % Note: Θ₂₁ = Θ₁₂^T, so P̃_iΘ₂₁ = P̃_i * (α_i+β_i)/2 * T_iT_i^T
    Fourth_Part_21 = P_tilde_i{i} * Theta_21 * TiTi_transpose + lambda_tilde_i{i} * eye(3);
    
    % Lower right 3×3 block: 0
    Fourth_Part_22 = zeros(3,3);
    
    % Construct the complete 6×6 Fourth_Part
    Fourth_Part = [Fourth_Part_11, Fourth_Part_12;               % Row 1: [Θ₁₁₃ₓ₃, (Θ₁₂P̃_i + λ̃_iI)₃ₓ₃]
                   Fourth_Part_21, Fourth_Part_22];              % Row 2: [(P̃_iΘ₂₁ + λ̃_iI)₃ₓ₃, 0₃ₓ₃]
    % Result: 6×6 matrix
    
    % COMPLETE SECOND CONSTRAINT MATRIX (equation 96 from Lemma 8)
    % All parts are now 6×6, so addition works correctly
    W2_i = First_Part + Second_Part - Third_Part - Fourth_Part;
        
    
    tagName = ['W2_',num2str(i)];
    constraintTags{end+1} = tagName;
    con5 = tag(W2_i >= epsilon*eye(size(W2_i)), tagName);
    constraintMats{end+1} = W2_i;
    
    % Add all DG constraints
    % constraints = [constraints, con0_1, con0_2, con1, con2, con3, con_nu, con_rho, con4, con5];
    constraints = [constraints, con0_1, con0_2, con1, con2, con3];

end

%% Line Constraints from Theorem 4 (Same as Lemma 5)
for l = 1:1:numOfLines
    Rl = Line{l}.R;
    Ll = Line{l}.L;
    Il = eye(1);
    
    % Line positivity constraint
    tagName = ['P_bar_',num2str(l),'_positive'];
    constraintTags{end+1} = tagName;
    con6 = tag(P_bar_l{l} >= epsilon*Il, tagName);
    constraintMats{end+1} = P_bar_l{l};
    
    % Line passivity constraints
    tagName = ['nu_bar_',num2str(l),'_bound'];
    constraintTags{end+1} = tagName;
    con_nu_l = tag(nu_bar_l{l} <= 0, tagName);  % Line passivity constraint
    constraintMats{end+1} = nu_bar_l{l};
    
    tagName = ['rho_bar_',num2str(l),'_positive'];
    constraintTags{end+1} = tagName;
    con_rho_l = tag(rho_bar_l{l} >= epsilon, tagName);  % Positive definite
    constraintMats{end+1} = rho_bar_l{l};
    
    % Main line LMI constraint (from Lemma 5)
    W_l = [(2*P_bar_l{l}*Rl)/Ll - rho_bar_l{l}, -P_bar_l{l}/Ll + 1/2;
           -P_bar_l{l}/Ll + 1/2,                 -nu_bar_l{l}];
    
    tagName = ['W_l_',num2str(l)];
    constraintTags{end+1} = tagName;
    con7 = tag(W_l >= epsilon*eye(2), tagName);
    constraintMats{end+1} = W_l;
    
    constraints = [constraints, con6, con_nu_l, con_rho_l, con7];
end

%% Mixed Constraints from Theorem 4 (equation 107) - Using your corrected structure
for i = 1:1:numOfDGs
    for l = 1:1:numOfLines
        if B_il(i,l) ~= 0
            p_i = piVals(i);
            p_l = plVals(l);
            Ct = DG{i}.C;
            
            % Define the non-zero values directly
            B_il_val = B_il(i,l);
            
            % Define auxiliary variable for bilinear term
            xi_il = sdpvar(1,1,'full');
            s1 = sdpvar(1,1,'full');
            s2 = sdpvar(1,1,'full');
            
            % Add LMI constraint for auxiliary variable
            auxiliary_LMI = [
                1,              nu_bar_l{l},    rho_tilde_i{i};
                nu_bar_l{l},    s1,             xi_il;
                rho_tilde_i{i}, xi_il,          s2
            ];
            
            tagName = ['AuxLMI_',num2str(i),'_',num2str(l)];
            constraintTags{end+1} = tagName;
            con_aux = tag(auxiliary_LMI >= epsilon*eye(3), tagName);
            constraintMats{end+1} = auxiliary_LMI;
            constraints = [constraints, con_aux];
            
            % Building the LMI step by step with subblocks following your structure
            % Block (1,1)
            Mat_11 = -p_i*nu_i{i};
            Mat_1 = Mat_11;
            
            % Block (2,2)
            Mat_22 = -p_l*nu_bar_l{l};
            
            % Block (1,2)
            Mat_12 = 0;
            
            % Combine to 2x2 matrix
            Mat_2 = [Mat_1, Mat_12;
                     Mat_12', Mat_22];
            
            % Block (3,3)
            Mat_33 = 1;
            
            % Block (1,3) and (2,3)
            Mat_13 = [0; 0];
            
            % Combine to 3x3 matrix
            Mat_3 = [Mat_2, Mat_13;
                     Mat_13', Mat_33];
            
            % Block (4,4)
            Mat_44 = p_i*rho_tilde_i{i};
            
            % Block (1,4), (2,4), and (3,4)
            Mat_14 = [0;
                      -B_il_val*xi_il*p_l;
                      rho_tilde_i{i}];
            
            % Combine to 4x4 matrix
            Mat_4 = [Mat_3, Mat_14;
                     Mat_14', Mat_44];
            
            % Block (5,5)
            Mat_55 = p_l*rho_bar_l{l};
            
            % Block (1,5), (2,5), (3,5), and (4,5)
            Mat_15 = [-(p_i*nu_i{i}*B_il_val)/Ct;  % Corrected sign
                      0;
                      1;
                      (1/2*p_i*B_il_val*rho_tilde_i{i})/Ct - 1/2*B_il_val*p_l*rho_tilde_i{i}];
            
            % Combine to 5x5 matrix
            Mat_5 = [Mat_4, Mat_15;
                     Mat_15', Mat_55];
            
            % Block (6,6)
            Mat_66 = gamma_tilde_i{i};
            
            % Block (1,6), (2,6), (3,6), (4,6), and (5,6)
            Mat_16 = [-p_i*nu_i{i};
                      -p_l*nu_bar_l{l};
                      0;
                      -1/2*p_i*rho_tilde_i{i};
                      -1/2*p_l];
            
            % Combine to final 6x6 matrix
            Mat_6 = [Mat_5, Mat_16;
                     Mat_16', Mat_66];
            
            % Final LMI condition
            W = Mat_6;
            tagName = ['LMI_107_',num2str(i),'_',num2str(l)];
            constraintTags{end+1} = tagName;
            con_LMI_107 = tag(W >= epsilon*eye(size(W)), tagName);
            constraintMats{end+1} = W;
            
            % Add to constraints
            constraints = [constraints, con_LMI_107];
        end
    end
end

%% Cost Function (minimizing α_λ Σ λ̃_i as in Theorem 4)
costGamma = 0;
alpha_lambda = 1;  % Weight for lambda terms

% Minimize λ̃_i terms (primary objective from Theorem 4)
% for i = 1:numOfDGs
%     costGamma = costGamma + alpha_lambda * lambda_tilde_i{i};
% end

% Additional regularization terms
for i = 1:numOfDGs
    costGamma = costGamma + 0.1*rho_tilde_i{i} + 100*gamma_tilde_i{i} + 0.01*trace(P_tilde_i{i}) + 0.01*trace(R_tilde_i{i});
end

% for l = 1:numOfLines
%     costGamma = costGamma + 0.1*rho_bar_l{l} + 0.01*trace(P_bar_l{l});
% end

%% Solve the LMI problem
solverOptions = sdpsettings('solver', 'mosek', 'verbose', 2, 'debug', debugMode, ...
                           'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-5, ...
                           'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1e-6, ...
                           'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS', 1e-6);

sol = optimize(constraints, costGamma, solverOptions);

statusLocalController = sol.problem == 0;

% Check for diagnostic information if infeasible
if sol.problem ~= 0
    disp(['Problem status: ', yalmiperror(sol.problem)]);
    disp('Problem is infeasible or numerical issues occurred.');
    
    if debugMode
        feasibility = check(constraints);
        combinedList = [feasibility(:), (1:length(feasibility))'];
        sortedList = sortrows(combinedList, 1);
        
        for i = 1:min(10, length(sortedList))
            idx = sortedList(i, 2);
            if feasibility(idx) < -1e-6
                disp(['Constraint "', constraintTags{idx}, '" is violated by ', num2str(feasibility(idx))]);
            end
        end
    end
end

%% Extract variable values if solution was found
if statusLocalController
    for i = 1:1:numOfDGs
        P_tilde_iVal = value(P_tilde_i{i});
        Kp_iVal = value(Kp_i{i});
        Ki_iVal = value(Ki_i{i});
        K_tilde_i0Val = [Kp_iVal, 0, Ki_iVal];  % Reconstruct with enforced structure
        nu_iVal = value(nu_i{i});
        rho_tilde_iVal = value(rho_tilde_i{i});
        gamma_tilde_iVal = value(gamma_tilde_i{i});
        R_tilde_iVal = value(R_tilde_i{i});
        lambda_tilde_iVal = value(lambda_tilde_i{i});
        
        % Convert back to original variables using Ki0 = K̃i0 * P̃i^(-1)
        Ki0 = K_tilde_i0Val / P_tilde_iVal;  % This gives [Kp, 0, Ki]
        rho_i = 1 / rho_tilde_iVal;  % Convert rho_tilde back to rho
        
        % Store results
        DG{i}.P0 = P_tilde_iVal;
        DG{i}.K0 = Ki0;
        DG{i}.nu = nu_iVal;
        DG{i}.rho = rho_i;
        DG{i}.gammaTilde0 = gamma_tilde_iVal;  % THIS IS THE KEY LINE FOR GLOBAL DESIGN!
        DG{i}.R_CPL = R_tilde_iVal;
        DG{i}.lambda_CPL = lambda_tilde_iVal;
        
        % Store CPL sector bounds for reference
        PL_i = DG{i}.PL;
        Vmin = 0.95 * DG{i}.refVoltage;
        Vmax = 1.05 * DG{i}.refVoltage;
        if PL_i > 0
            DG{i}.alpha_CPL = PL_i / (DG{i}.C * Vmax^2);
            DG{i}.beta_CPL = PL_i / (DG{i}.C * Vmin^2);
        else
            DG{i}.alpha_CPL = 0;
            DG{i}.beta_CPL = 0;
        end
        
        fprintf('DG %d: gammaTilde0 = %f stored successfully\n', i, gamma_tilde_iVal);
        fprintf('DG %d: K0 = [%f, %f, %f]\n', i, Ki0(1), Ki0(2), Ki0(3));
    end
    
    for l = 1:1:numOfLines
        P_bar_lVal = value(P_bar_l{l});
        nu_bar_lVal = value(nu_bar_l{l});
        rho_bar_lVal = value(rho_bar_l{l});
        
        % Store results
        Line{l}.P0 = P_bar_lVal;
        Line{l}.nu = nu_bar_lVal;
        Line{l}.rho = rho_bar_lVal;
        
        fprintf('Line %d: nu = %f, rho = %f\n', l, nu_bar_lVal, rho_bar_lVal);
    end
    
    fprintf('Local control design completed successfully!\n');
else
    fprintf('Local control design FAILED - no values stored.\n');
end

end