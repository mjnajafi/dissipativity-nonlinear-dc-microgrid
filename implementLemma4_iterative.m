function [Vr, Is, statusLemma4] = implementLemma4_iterative(DG, Line, B_il, maxIter)
    if nargin < 4
        maxIter = 10;
    end
    
    numDGs = length(DG);
    numLines = length(Line);
    
    % Extract system parameters
    Rt = zeros(numDGs, numDGs);
    YL = zeros(numDGs, numDGs);
    Pn = zeros(numDGs, 1);
    IL_bar = zeros(numDGs, 1);
    PL = zeros(numDGs, 1);  % Constant power loads
    V_nominal = zeros(numDGs, 1);
    
    for i = 1:numDGs
        Rt(i,i) = DG{i}.R;
        YL(i,i) = DG{i}.Y;
        V_nominal(i) = DG{i}.refVoltage;
        Pn(i) = DG{i}.powerRating / DG{i}.refVoltage;
        IL_bar(i) = DG{i}.IL;
        
        % Extract constant power load
        if isfield(DG{i}, 'PL')
            PL(i) = DG{i}.PL;
        else
            PL(i) = 0;  % No CPL
        end
    end
    
    % Build R matrix for lines
    R = zeros(numLines, numLines);
    for l = 1:numLines
        R(l,l) = Line{l}.R;
    end
    
    % Voltage bounds
    voltage_deviation_percent = 10;  % 10% deviation
    Vmin = (1 - voltage_deviation_percent/100) * V_nominal;
    Vmax = (1 + voltage_deviation_percent/100) * V_nominal;
    
    % Initialize with nominal voltages
    Vr_prev = V_nominal;
    tolerance = 1e-6;
    
    fprintf('\n==== Implementing Lemma 4 with Iterative Linearization ====\n');
    
    for iter = 1:maxIter
        fprintf('Iteration %d:\n', iter);
        
        % Linearize diag(Vr)^(-1)*PL around Vr_prev
        % f(Vr) ≈ f(Vr_prev) + ∇f(Vr_prev) * (Vr - Vr_prev)
        % where f(Vr) = diag(Vr)^(-1)*PL
        
        % Current linearization point
        f_prev = PL ./ Vr_prev;  % Element-wise division
        
        % Gradient: ∂/∂Vr_i [PL_i/Vr_i] = -PL_i/Vr_i^2
        grad_diag = -PL ./ (Vr_prev.^2);
        
        % Linear approximation: f(Vr) ≈ f_prev + grad_diag .* (Vr - Vr_prev)
        % This gives us: PL./Vr ≈ f_prev + grad_diag .* (Vr - Vr_prev)
        
        % Define optimization variables
        Vr = sdpvar(numDGs, 1);
        Is = sdpvar(1, 1);
        
        % Linear approximation of CPL current
        CPL_current_linear = f_prev + grad_diag .* (Vr - Vr_prev);
        
        % Create constraints
        constraints = [];
        
        % Main equilibrium constraint with linearized CPL term
        % ItE = (B*R^(-1)*B' + YL)*Vr + IL_bar + CPL_current_linear = Pn*Is
        equilibrium_current = (B_il*inv(R)*B_il' + YL)*Vr + IL_bar + CPL_current_linear;
        constraints = [constraints, equilibrium_current == Pn*Is];
        
        % Bounds
        constraints = [constraints, Vmin <= Vr <= Vmax];
        constraints = [constraints, 0 <= Is <= 1];
        
        % Objective: minimize deviation from nominal + current sharing coefficient
        alphaV = 1;
        alphaI = 0.1;
        objective = alphaV * norm(Vr - V_nominal, 2)^2 + alphaI * Is^2;
        
        % Solve
        options = sdpsettings('solver', 'mosek', 'verbose', 0);
        sol = optimize(constraints, objective, options);
        
        if sol.problem ~= 0
            fprintf('  Optimization failed at iteration %d\n', iter);
            statusLemma4 = false;
            return;
        end
        
        % Extract solution
        Vr_new = value(Vr);
        Is_new = value(Is);
        
        % Check convergence
        voltage_change = norm(Vr_new - Vr_prev);
        fprintf('  Voltage change: %.6f\n', voltage_change);
        
        if voltage_change < tolerance
            fprintf('  Converged!\n');
            Vr = Vr_new;
            Is = Is_new;
            statusLemma4 = true;
            break;
        end
        
        % Update for next iteration
        Vr_prev = Vr_new;
        
        if iter == maxIter
            fprintf('  Maximum iterations reached\n');
            Vr = Vr_new;
            Is = Is_new;
            statusLemma4 = false;  % Not converged
        end
    end
    
    if statusLemma4
        % Verify the solution
        fprintf('\n==== Final Solution Verification ====\n');
        fprintf('Current sharing coefficient Is: %.4f\n', Is);
        
        % Compute actual equilibrium currents (with nonlinear CPL)
        actual_CPL_current = PL ./ Vr;  % Actual nonlinear term
        actual_currents = (B_il*inv(R)*B_il' + YL)*Vr + IL_bar + actual_CPL_current;
        expected_currents = Pn * Is;
        
        fprintf('\nCurrent sharing verification:\n');
        for i = 1:numDGs
            fprintf('  DG%d: Actual=%.3f A, Expected=%.3f A, Error=%.6f\n', ...
                i, actual_currents(i), expected_currents(i), ...
                abs(actual_currents(i) - expected_currents(i)));
        end
        
        % Check equilibrium error
        equilibrium_error = norm(actual_currents - expected_currents);
        fprintf('Total equilibrium error: %.6f\n', equilibrium_error);
        
        if equilibrium_error > 1e-3
            fprintf('Warning: Large equilibrium error detected!\n');
        end
    end
end