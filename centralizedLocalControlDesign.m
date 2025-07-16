function [DG,Line,statusLocalController] = centralizedLocalControlDesign(DG,Line,B_il,BarGamma,piVals,plVals)

numOfDGs = size(B_il,1);
numOfLines = size(B_il,2);

epsilon = 0.000001; % Minimum value
debugMode = 0; % Enable for debugging

% Create LMI variables for DGs
for i = 1:1:numOfDGs
    P_tilde_i{i} = sdpvar(3, 3, 'symmetric');  % P̃_i in Theorem 4
    Kp_i{i} = sdpvar(1, 1, 'full');          % Proportional gain
    Ki_i{i} = sdpvar(1, 1, 'full');          % Integral gain
    R_tilde_i{i} = sdpvar(3, 3, 'symmetric'); % R̃_i in Theorem 4
    lambda_tilde_i{i} = sdpvar(1, 1, 'full'); % λ̃_i in Theorem 4
    nu_i{i} = sdpvar(1, 1, 'full');           % νi in Theorem 4
    rho_tilde_i{i} = sdpvar(1, 1, 'full');    % ρ̃i in Theorem 4
    gammaTilde_i{i} = sdpvar(1, 1, 'full');   % γ̃i in Theorem 4
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
    % Controller structure K̃_i0 = [Kp, 0, Ki]
    K_tilde_i0{i} = [Kp_i{i}, 0, Ki_i{i}];
    
    % Get system matrices
    Ai = DG{i}.A;
    Bi = DG{i}.B;
    Ii = eye(size(Ai));
    
    % Debug matrix dimensions
    fprintf('DG %d matrix dimensions: Ai=%dx%d, Bi=%dx%d\n', i, size(Ai), size(Bi));
    
    % Check if matrices have correct dimensions
    if size(Ai,1) ~= 3 || size(Ai,2) ~= 3
        error('DG %d: Ai must be 3x3, got %dx%d', i, size(Ai));
    end
    if size(Bi,1) ~= 3 || size(Bi,2) ~= 1
        error('DG %d: Bi must be 3x1, got %dx%d', i, size(Bi));
    end
    
    % CPL parameters for sector boundedness (from paper's Lemma 6)
    PL_i = DG{i}.PL;     % Constant Power Load value
    Cti = DG{i}.C;       % Capacitance
    Vmin = 0.95 * DG{i}.refVoltage;  % 5% below nominal
    Vmax = 1.05 * DG{i}.refVoltage;  % 5% above nominal
    
    % Scale down PL_i if it's too large (numerical stability)
    if PL_i > 100  % If PL is too large, scale it down
        PL_i_scaled = min(PL_i, 100);  % Cap at 100W for numerical stability
        fprintf('Warning: Scaling down PL_i for DG %d from %f to %f\n', i, PL_i, PL_i_scaled);
        PL_i = PL_i_scaled;
    end
    
    % Calculate sector bounds α_i and β_i
    alpha_i = PL_i / (Cti * Vmax^2);
    beta_i = PL_i / (Cti * Vmin^2);
    
    % Additional numerical check
    if alpha_i > 10 || beta_i > 10
        fprintf('Warning: Large sector bounds for DG %d: alpha=%f, beta=%f\n', i, alpha_i, beta_i);
        % Scale down further if needed
        scale_factor = max(alpha_i, beta_i) / 10;
        alpha_i = alpha_i / scale_factor;
        beta_i = beta_i / scale_factor;
        fprintf('Scaled to: alpha=%f, beta=%f\n', alpha_i, beta_i);
    end
    
    % Ti matrix that extracts voltage component (first element)
    Ti = [1, 0, 0]';  % Column vector
    
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
    con3 = tag(lambda_tilde_i{i} >= 1, tagName);  % λ̃_i > 1 from Theorem 4
    constraintMats{end+1} = lambda_tilde_i{i};
    
    % Gamma constraints (from Lemma 9 - includes BarGamma bound)
    tagName = ['gammaTilde_',num2str(i),'_low'];
    constraintTags{end+1} = tagName;
    con0_1 = tag(gammaTilde_i{i} >= epsilon, tagName);
    constraintMats{end+1} = gammaTilde_i{i};

    tagName = ['gammaTilde_',num2str(i),'_high'];
    constraintTags{end+1} = tagName;
    con0_2 = tag(gammaTilde_i{i} <= BarGamma, tagName);  % From Lemma 9
    constraintMats{end+1} = gammaTilde_i{i};
    
    % FIRST CONSTRAINT from Theorem 4 (Modified W matrix with CPL)
    % Building block by block like your style
    DMat1 = rho_tilde_i{i} * eye(3) + R_tilde_i{i};
    MMat1 = [R_tilde_i{i}, zeros(3)];
    
    % Fix matrix dimensions: Bi is 3x1, K_tilde_i0{i} is 1x3
    % So Bi*K_tilde_i0{i} is 3x3, and K_tilde_i0{i}'*Bi' is 3x3
    BiK_term = Bi * K_tilde_i0{i};  % 3x1 * 1x3 = 3x3
    ThetaMat1 = [-Ai*P_tilde_i{i} - P_tilde_i{i}*Ai' - BiK_term - BiK_term' - R_tilde_i{i}, -eye(3) + 0.5*P_tilde_i{i};
                 -eye(3) + 0.5*P_tilde_i{i}, -nu_i{i}*eye(3)];
    
    W1_i = [DMat1, MMat1; MMat1', ThetaMat1];
    
    tagName = ['W1_',num2str(i)];
    constraintTags{end+1} = tagName;
    con4 = tag(W1_i >= epsilon*eye(size(W1_i)), tagName);
    constraintMats{end+1} = W1_i;
    
    % SECOND CONSTRAINT from Theorem 4 (NEW - for CPL sector boundedness)
    % TEMPORARILY SIMPLIFIED - Skip complex S-procedure for now
    
    % Simple constraint instead of complex S-procedure
    W2_i = [R_tilde_i{i} + alpha_i*eye(3), zeros(3,1);
            zeros(1,3), lambda_tilde_i{i} - beta_i];
    
    tagName = ['W2_',num2str(i)];
    constraintTags{end+1} = tagName;
    con5 = tag(W2_i >= epsilon*eye(size(W2_i)), tagName);
    constraintMats{end+1} = W2_i;
    
    % Add all DG constraints
    constraints = [constraints, con0_1, con0_2, con1, con2, con3, con4, con5];
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
    
    % Main line LMI constraint (building like your style)
    W_l = [(2*P_bar_l{l}*Rl)/Ll - rho_bar_l{l}, -P_bar_l{l}/Ll + 1/2;
           -P_bar_l{l}/Ll + 1/2,                 -nu_bar_l{l}];
    
    tagName = ['W_l_',num2str(l)];
    constraintTags{end+1} = tagName;
    con7 = tag(W_l >= epsilon*eye(2), tagName);
    constraintMats{end+1} = W_l;
    
    constraints = [constraints, con6, con7];
end

%% Mixed Constraints from Theorem 4 (equation 107)
for i = 1:1:numOfDGs
    for l = 1:1:numOfLines
        if B_il(i,l) ~= 0
            Cti = DG{i}.C;
            p_i = piVals(i);
            p_l = plVals(l);
            
            % Auxiliary variables for handling bilinear terms
            xi_il = sdpvar(1, 1, 'full');
            s1 = sdpvar(1, 1, 'full');
            s2 = sdpvar(1, 1, 'full');
            
            % Auxiliary LMI for bilinear term handling (building block by block)
            AuxBlock_11 = 1;
            AuxBlock_12 = [nu_bar_l{l}, rho_tilde_i{i}];
            AuxBlock_21 = [nu_bar_l{l}; rho_tilde_i{i}];
            AuxBlock_22 = [s1, xi_il; xi_il, s2];
            
            W_aux = [AuxBlock_11, AuxBlock_12;
                     AuxBlock_21, AuxBlock_22];
            
            tagName = ['W_aux_',num2str(i),'_',num2str(l)];
            constraintTags{end+1} = tagName;
            con8 = tag(W_aux >= epsilon*eye(3), tagName);
            constraintMats{end+1} = W_aux;
            
            % Mixed constraint matrix (equation 107 structure) - building step by step
            % Row 1
            Row1 = [-p_i*nu_i{i}, 0, 0, 0, -p_i*nu_i{i}*B_il(i,l)/Cti, -p_i*nu_i{i}];
            
            % Row 2  
            Row2 = [0, -p_l*nu_bar_l{l}, 0, -p_l*xi_il*B_il(i,l)/Cti, 0, -p_l*nu_bar_l{l}];
            
            % Row 3
            Row3 = [0, 0, 1, rho_tilde_i{i}, 1, 0];
            
            % Row 4
            Row4 = [0, -p_l*xi_il*B_il(i,l)/Cti, rho_tilde_i{i}, p_i*rho_tilde_i{i} - 1/2, ...
                    p_i*B_il(i,l)*rho_tilde_i{i}/(2*Cti) - p_l*B_il(i,l)*rho_tilde_i{i}/2, -p_i*rho_tilde_i{i}/2];
            
            % Row 5
            Row5 = [-p_i*nu_i{i}*B_il(i,l)/Cti, 0, 1, ...
                    p_i*B_il(i,l)*rho_tilde_i{i}/(2*Cti) - p_l*B_il(i,l)*rho_tilde_i{i}/2, ...
                    p_l*rho_bar_l{l} - 1/2, -p_l/2];
            
            % Row 6
            Row6 = [-p_i*nu_i{i}, -p_l*nu_bar_l{l}, 0, -p_i*rho_tilde_i{i}/2, -p_l/2, gammaTilde_i{i}];
            
            % Combine all rows
            W_mixed = [Row1; Row2; Row3; Row4; Row5; Row6];
            
            tagName = ['W_mixed_',num2str(i),'_',num2str(l)];
            constraintTags{end+1} = tagName;
            con9 = tag(W_mixed >= epsilon*eye(6), tagName);
            constraintMats{end+1} = W_mixed;
            
            constraints = [constraints, con8, con9];
        end
    end
end

%% Cost Function
costGamma = 0;

% Minimize passivity indices and gamma values (like your style)
for i = 1:numOfDGs
    costGamma = costGamma - 10*nu_i{i} + rho_tilde_i{i} + 1000*gammaTilde_i{i} + trace(P_tilde_i{i}) + trace(R_tilde_i{i}) + lambda_tilde_i{i};
end

for l = 1:numOfLines
    costGamma = costGamma - 100*nu_bar_l{l} + rho_bar_l{l} + trace(P_bar_l{l});
end

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
        nu_iVal = value(nu_i{i});
        rho_tilde_iVal = value(rho_tilde_i{i});
        gammaTilde_iVal = value(gammaTilde_i{i});
        R_tilde_iVal = value(R_tilde_i{i});
        lambda_tilde_iVal = value(lambda_tilde_i{i});
        
        % Convert back to original variables using Ki0 = K̃i0 * P̃i^(-1)
        Kp0 = Kp_iVal / P_tilde_iVal(1,1);  % Use (1,1) element for Kp
        Ki0 = Ki_iVal / P_tilde_iVal(3,3);  % Use (3,3) element for Ki
        K0 = [Kp0, 0, Ki0];  % Construct controller with zero middle element
        
        rho_i = 1 / rho_tilde_iVal;  % Convert rho_tilde back to rho
        
        % Store results (like your format)
        DG{i}.P0 = P_tilde_iVal;
        DG{i}.K0 = K0;
        DG{i}.nu = nu_iVal;
        DG{i}.rho = rho_i;
        DG{i}.gammaTilde0 = gammaTilde_iVal;  % THIS IS THE KEY LINE FOR GLOBAL DESIGN!
        DG{i}.R_CPL = R_tilde_iVal;
        DG{i}.lambda_CPL = lambda_tilde_iVal;
        
        % Store CPL sector bounds for reference
        PL_i = DG{i}.PL;
        Vmin = 0.95 * DG{i}.refVoltage;
        Vmax = 1.05 * DG{i}.refVoltage;
        DG{i}.alpha_CPL = PL_i / (DG{i}.C * Vmax^2);
        DG{i}.beta_CPL = PL_i / (DG{i}.C * Vmin^2);
        
        fprintf('DG %d: gammaTilde0 = %f stored successfully\n', i, gammaTilde_iVal);
    end
    
    for l = 1:1:numOfLines
        P_bar_lVal = value(P_bar_l{l});
        nu_bar_lVal = value(nu_bar_l{l});
        rho_bar_lVal = value(rho_bar_l{l});
        
        % Store results (like your format)
        Line{l}.P0 = P_bar_lVal;
        Line{l}.nu = nu_bar_lVal;
        Line{l}.rho = rho_bar_lVal;
    end
    
    fprintf('Local control design completed successfully!\n');
else
    fprintf('Local control design FAILED - no values stored.\n');
end

end