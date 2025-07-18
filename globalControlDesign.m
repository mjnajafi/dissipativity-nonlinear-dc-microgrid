% This function constructs the control matrices C, BarC, and H, and 
% formulates the optimization problem with appropriate constraints and 
% objective functions to ensure the stability and performance of the 
% microgrid under a given communication topology defined by the adjacency 
% matrix A_ij. The function also computes the optimal gain matrices K for 
% each pair of DGs based on the optimization results.

function [DG,Line,statusGlobalController,gammaTildeVal,K,C,BarC,H,P_iVal,P_lVal] = globalControlDesign(DG,Line,A_ij,B_il,BarGamma,isSoft)

% Add verification step
% verificationStatus = verifyControlConditions(DG, Line, B_il);
% if ~verificationStatus
%     fprintf('WARNING: Control design conditions not satisfied.\n');
%     fprintf('Global controller design may fail.\n\n');
% end


numOfDGs = size(B_il,1);
numOfLines = size(B_il,2);
epsilon = 0.000001; % Minimum value
debugMode = 0;

%% Creating C , BarC , and H Matrices

% Create C Matrix
C = zeros(numOfLines, numOfDGs * 3);

for l = 1:numOfLines
    for i = 1:numOfDGs
        C(l, (i-1) * 3 + 1) = B_il(i, l);
    end
end



% Create BarC Matrix
BarC = zeros(numOfDGs * 3, numOfLines);

for i = 1:numOfDGs
    for l = 1:numOfLines
        Ct = DG{i}.C;
        BarC((i-1)*3 + 1, l) = -B_il(i, l)/Ct;
    end
end


% Create H Matrix
H_i = [0, 1 , 1];
H = [];
for i = 1:numOfDGs
    H = blkdiag(H, H_i);
end

% %this is debugging
% fprintf('Debugging DG structure:\n');
% for i = 1:numOfDGs
%     fprintf('DG{%d} fields: ', i);
%     disp(fieldnames(DG{i}));
% 
%     if isfield(DG{i}, 'gammaTilde0')
%         fprintf('DG{%d}.gammaTilde0 = %f\n', i, DG{i}.gammaTilde0);
%     else
%         error('ERROR: DG{%d}.gammaTilde0 does not exist!', i);
%     end
% end



% GammaMat
GammaMat = [];
for i = 1:numOfDGs
    GammaMat = blkdiag(GammaMat, DG{i}.gammaTilde0);
end


%% Creating the adjacency matrix, null matrix, and cost matrix
A = A_ij; % Adjacency matrix of the DG-DG communication topology - Physical Topology
% A = ones(numOfDGs) - eye(numOfDGs); % all-to-all communication topology
% A = zeros(numOfDGs); % all-to-all communication topology



adjMatBlock = cell(numOfDGs, numOfDGs);
nullMatBlock = cell(numOfDGs, numOfDGs);
costMatBlock = cell(numOfDGs, numOfDGs);

for i = 1:numOfDGs
    for j = 1:numOfDGs
        if i ~= j  % Only for non-diagonal terms
            if A(j, i) == 1
                % Structure for direct difference gains
                adjMatBlock{i, j} = [0, 0, 0; 0, 1, 0; 0, 0, 0];  % This will be the gain for (Ij - Ii)
                nullMatBlock{i, j} = [1, 1, 1; 1, 0, 1; 1, 1, 1];
                % Cost based on physical distance
                dist_ij = norm(DG{i}.coordinates - DG{j}.coordinates);
                costMatBlock{i, j} = dist_ij^2 * [0, 0, 0; 0, 1, 0; 0, 0, 0];
            else
                % Non-connected DGs
                adjMatBlock{i, j} = zeros(3);
                nullMatBlock{i, j} = [1, 1, 1; 1, 0, 1; 1, 1, 1];
                costMatBlock{i, j} = zeros(3);
            end
        else
            % Diagonal terms (i == j) - all zero as we don't need self-gains
            adjMatBlock{i, j} = zeros(3);
            nullMatBlock{i, j} = ones(3);  % Force diagonal blocks to be zero
            costMatBlock{i, j} = zeros(3);
        end
    end
end

%%%% Comment: Use the following variable.
adjMatBlock = cell2mat(adjMatBlock);  %%% This is only for enforcing hard garaph constraints
nullMatBlock = cell2mat(nullMatBlock); %%% This is needed in any case as it is a structural constraint
costMatBlock = cell2mat(costMatBlock);   %%% Note that this matrix contains the c_ij information required for the objective function (46a)
%% Variables corresponding to DGs like 
I = eye(3*numOfDGs);
I_N = eye(numOfDGs);
gammaTilde = sdpvar(1, 1,'full');
GammaTilde = gammaTilde*I;

Q = sdpvar(3*numOfDGs, 3*numOfDGs, 'full'); 

P_i = sdpvar(numOfDGs, numOfDGs, 'diagonal');

P_l = sdpvar(numOfLines, numOfLines, 'diagonal');

I_n = eye(3);


X_p_11 = [];
X_p_12 = [];
X_12 = [];
X_p_22 = [];
for i = 1:1:numOfDGs
          
    nu_i = DG{i}.nu;
    rho_i = DG{i}.rho;

    % For DGs, according to equation (33), numOfInouts = numOfOutputs = 3
    X_i_11 = -nu_i*I_n;    %inputs x inputs
    X_i_12 = 0.5*I_n;      %inputs x outputs
    X_i_22 = -rho_i*I_n;   %outputs x outputs

    X_p_11 = blkdiag(X_p_11, P_i(i,i)*X_i_11);
    X_p_12 = blkdiag(X_p_12, P_i(i,i)*X_i_12);
    X_p_22 = blkdiag(X_p_22, P_i(i,i)*X_i_22);
    X_12 = blkdiag(X_12, (-1 / (2 * nu_i)) * I_n);
end
X_p_21 = X_p_12';
X_21 = X_12';

I_bar = eye(1);

BarX_Barp_11 = [];
BarX_Barp_12 = [];
BarX_12 = [];
BarX_Barp_22 = [];

for l = 1:1:numOfLines
    
    nu_l = Line{l}.nu;
    rho_l = Line{l}.rho;

    % For Lines, according to equation (31), numOfInputs = numOfOutputs = 1
    BarX_l_11 = -nu_l*I_bar;    %inputs x inputs
    BarX_l_12 = 0.5*I_bar;      %inputs x outputs
    BarX_l_22 = -rho_l*I_bar;   %outputs x outputs

    BarX_Barp_11 = blkdiag(BarX_Barp_11, P_l(l,l)*BarX_l_11);
    BarX_Barp_12 = blkdiag(BarX_Barp_12, P_l(l,l)*BarX_l_12);
    BarX_Barp_22 = blkdiag(BarX_Barp_22, P_l(l,l)*BarX_l_22);
    
    BarX_12 = blkdiag(BarX_12, (-1 / (2 * nu_l)) * I_bar);
    
end
BarX_Barp_21 = BarX_Barp_12';
BarX_21 = BarX_12';



%% Constraints 
constraints = [];
constraintTags = {}; % Cell array to hold tags
constraintMats = {}; % Cell array to hold matrices

constraints = [constraints,
    norm(Q, 'fro') <= 1e4,
    norm(P_i, 'fro') <= 1e4,
    norm(P_l, 'fro') <= 1e4,
    norm(GammaTilde, 'fro') <= 1e4
];

% Constraints in (46b)
for i = 1:numOfDGs
    tagName = ['P_',num2str(i),num2str(i)];
    constraintTags{end+1} = tagName;
    con1 = tag(P_i(i,i) >= epsilon, tagName);
    constraintMats{end+1} = P_i(i,i);

    constraints = [constraints, con1];
end

% Constraints in (46c)
for l = 1:numOfLines
    tagName = ['PBar_',num2str(l),num2str(l)];
    constraintTags{end+1} = tagName;
    con2 = tag(P_l(l,l) >= epsilon, tagName);
    constraintMats{end+1} = P_l(l,l);

    constraints = [constraints, con2];
end

% Constraints in (46d)
tagName = ['gammaTilde_low'];
constraintTags{end+1} = tagName;
con3_1 = tag(gammaTilde >= epsilon, tagName);
constraintMats{end+1} = gammaTilde;

tagName = ['gammaTilde_high'];
constraintTags{end+1} = tagName;
con3_2 = tag(gammaTilde <= BarGamma, tagName);
constraintMats{end+1} = gammaTilde;

constraints = [constraints, con3_1, con3_2];

% Constraint in (47)
O_n = zeros(numOfDGs, 3*numOfDGs);
O_bar = zeros(numOfLines);
O_N = zeros(3*numOfDGs, numOfLines);
O = zeros(numOfDGs, numOfLines);
O_3N = zeros(3*numOfDGs);

Mat_11 = X_p_11; 
Mat_1 = [Mat_11];
 

Mat_22 = BarX_Barp_11;
Mat_12 = O_N;
Mat_2 = [Mat_1, Mat_12;
        Mat_12', Mat_22];


Mat_33 = I_N;
Mat_13 = [  O_n';
            O'      ];
Mat_3 = [Mat_2, Mat_13;
        Mat_13', Mat_33];




Mat_44 = -Q' * X_12 - X_21 * Q - X_p_22;
Mat_14 = [  Q;  
            BarX_Barp_11 * C;  
            H                   ];
Mat_4 = [Mat_3, Mat_14;
        Mat_14', Mat_44];

Mat_4_Test1 = [ Mat_11, Q;
                Q', Mat_44];
Mat_4_Test2 = [ Mat_22, BarX_Barp_11 * C;
                (BarX_Barp_11 * C)', Mat_44];
Mat_4_Test3 = [ Mat_33, H;
                H', Mat_44];


Mat_55 =  -BarX_Barp_22;
Mat_15 = [  X_p_11 * BarC;
            O_bar;
            O;   
            -X_21 * X_p_11 * BarC - C' * BarX_Barp_11' * BarX_12];
Mat_5 = [Mat_4, Mat_15;
        Mat_15', Mat_55];

Mat_5_Test1 = [Mat_11, (X_p_11 * BarC);
               (X_p_11 * BarC)', Mat_55];
Mat_5_Test2 = [Mat_44, (-X_21 * X_p_11 * BarC - C' * BarX_Barp_11' * BarX_12);
               (-X_21 * X_p_11 * BarC - C' * BarX_Barp_11' * BarX_12)', Mat_55];

Mat_66 = GammaTilde;
Mat_16 = [  X_p_11;
            O_N';
            O_n;
            -X_21 * X_p_11;
            O_N'    ];
Mat_6 = [Mat_5, Mat_16;
        Mat_16', Mat_66];

Mat_6_Test1 = [Mat_11, X_p_11;
                X_p_11', Mat_66];
Mat_6_Test2 = [Mat_44, (-X_21 * X_p_11);
                (-X_21 * X_p_11)', Mat_66];


W = Mat_6;

slack_W = sdpvar(size(W,1), size(W,1), 'symmetric');
slack_mat44 = sdpvar(size(Mat_44,1), size(Mat_44,1), 'symmetric');
slack_mat55 = sdpvar(size(Mat_55,1), size(Mat_55,1), 'symmetric');
slack_test1 = sdpvar(size(Mat_6_Test1,1), size(Mat_6_Test1,1), 'symmetric');
slack_test2 = sdpvar(size(Mat_6_Test2,1), size(Mat_6_Test2,1), 'symmetric');

tagName = ['W'];
constraintTags{end+1} = tagName;
con4 = tag(W + slack_W >= epsilon*eye(size(W)), 'W');
% con4 = tag(W >= epsilon*eye(size(W)), tagName);
constraintMats{end+1} = W;
constraints = [constraints, con4];


constraints = [constraints,
    Mat_6_Test1 + slack_test1 >= epsilon*eye(size(Mat_6_Test1)),
    Mat_6_Test2 + slack_test2 >= epsilon*eye(size(Mat_6_Test2)),
    slack_test1 >= 0,
    slack_test2 >= 0,
    trace(slack_test1) <= 1e3,
    trace(slack_test2) <= 1e3,
    Mat_44 + slack_mat44 >= -1e4*eye(size(Mat_44)),
    Mat_55 + slack_mat55 >= -1e4*eye(size(Mat_55)),
    slack_W >= 0,
    slack_mat44 >= 0,
    slack_mat55 >= 0,
    trace(slack_W) <= 1e3,      % Limit the total slack
    trace(slack_mat44) <= 1e3,
    trace(slack_mat55) <= 1e3
];


% % Add these diagnostic checks after constructing Mat_6
% fprintf('Matrix Properties Check:\n');
% fprintf('Size of Mat_6: %dx%d\n', size(Mat_6));
% fprintf('Is Mat_6 symmetric? %d\n', norm(Mat_6 - Mat_6') < 1e-10);
% 
% % Check individual block properties
% fprintf('\nBlock Matrix Properties:\n');
% fprintf('Mat_44 eigenvalues range: [%.2e, %.2e]\n', min(eig(value(Mat_44))), max(eig(value(Mat_44))));
% fprintf('Mat_55 eigenvalues range: [%.2e, %.2e]\n', min(eig(value(Mat_55))), max(eig(value(Mat_55))));
% 



% % Temporary - 1:
% W_Temp = Mat_6_Test1;
% tagName = ['W_Temp_Mat_6_Test1'];
% constraintTags{end+1} = tagName;
% con4_Temp = tag(W_Temp >= epsilon*eye(size(W_Temp)), tagName);
% constraintMats{end+1} = W_Temp;
% constraints = [constraints, con4_Temp];
% 
% % Temporary - 2:
% W_Temp = Mat_6_Test2;
% tagName = ['W_Temp_Mat_6_Test2'];
% constraintTags{end+1} = tagName;
% con4_Temp = tag(W_Temp >= epsilon*eye(size(W_Temp)), tagName);
% constraintMats{end+1} = W_Temp;
% constraints = [constraints, con4_Temp];

% % Temporary - 3:
% W_Temp = Mat_5_Test1;
% tagName = ['W_Temp_Mat_5_Test1'];
% constraintTags{end+1} = tagName;
% con4_Temp = tag(W_Temp >= epsilon*eye(size(W_Temp)), tagName);
% constraintMats{end+1} = W_Temp;
% constraints = [constraints, con4_Temp];
% 
% % Temporary - 4:
% W_Temp = Mat_5_Test2;
% tagName = ['W_Temp_Mat_5_Test2'];
% constraintTags{end+1} = tagName;
% con4_Temp = tag(W_Temp >= epsilon*eye(size(W_Temp)), tagName);
% constraintMats{end+1} = W_Temp;
% constraints = [constraints, con4_Temp];

% 49g and h (Mat_5_Test1 & 2) conflicts 49e (Mat_4_test1 & 2)




% Structural constraints
tagName = ['Q_Structure'];
constraintTags{end+1} = tagName;
con5 = tag(Q.*(nullMatBlock==1) == O_3N, tagName);     % Structural limitations (due to the format of the control law)
constraintMats{end+1} = Q;

constraints = [constraints, con5];





% Hard Graph Constraints (forcing K_ij = K_ji = 0 if i and j are not communication neighbors)
if ~isSoft
    % Follow strictly the given communication tpology (by A_ij adjacency matrix)
    tagName = ['Q_Topology'];
    constraintTags{end+1} = tagName;
    con7 = tag(Q.*(adjMatBlock==0) == O_3N, tagName);      % Graph structure : hard constraint
    constraintMats{end+1} = Q;

    constraints = [constraints, con7]; % With the hard graph constraint con7
end



%% Objective Function
normType = 2;
% costFun0 = 1*norm(Q.*costMatBlock,normType);
costFun0 = sum(sum(Q.*costMatBlock)); %%% Play with this

% Additional constraint on costFun0
% Minimum Budget Constraints
con6 = costFun0 >= 0.01;           % Add the con6 constraint
constraints = [constraints, con6];  % Include con6 in the constraints

if isSoft
    alpha = 100;  % When isSoft is 1, emphasize communication cost
else
    alpha = 0;  % When isSoft is 0, exclude communication cost
end

beta = 0.1;       % Robustness - performance cost

% Total Cost Function
% costFun = alpha * costFun0 + beta * gammaTilde;
% Add slack penalties to the cost function
slackPenalty = 1000;  % Weight for slack terms
costFun = alpha * costFun0 + beta * gammaTilde + ...
    slackPenalty * (trace(slack_W) + trace(slack_mat44) + trace(slack_mat55) + trace(slack_test1) + trace(slack_test2));


%% Solve the LMI problem (47)

solverOptions = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 0);
% solverOptions = sdpsettings('solver', 'mosek', ...
%     'verbose', 1, ...
%     'debug', 1, ...
%     'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-8, ...
%     'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1e-8, ...
%     'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS', 1e-8, ...
%     'mosek.MSK_DPAR_PRESOLVE_TOL_X', 1e-8);

sol = optimize(constraints,costFun,solverOptions);

statusGlobalController = sol.problem == 0;   

%% Extract variable values

P_iVal = diag(value(P_i));
P_lVal = diag(value(P_l));

gammaTildeVal = value(gammaTilde)
costFun0Val = value(costFun0)

QVal = value(Q);
X_p_11Val = value(X_p_11);
KVal = X_p_11Val \ QVal;




% Load KVal elements into a cell structure K{i,j} (i.e., partitioning KVal
% into N\times N blocks)
% Obtaining K_ij blocks
KVal(nullMatBlock == 1) = 0;
maxNorm = 0;
for i = 1:1:numOfDGs
    for j = 1:1:numOfDGs
        K{i,j} = KVal(3*(i-1)+1:3*i , 3*(j-1)+1:3*j); % (i,j)-th (3 x 3) block
        normVal = max(max(abs(K{i,j})));
        if normVal>maxNorm 
            maxNorm = normVal;
        end
    end
end

% % Filter the K_ij values (weed out the ones with the smallest magnitudes)
% filtering out extremely small interconnections
% for i=1:1:numOfDGs
%     for j=1:1:numOfDGs
%         if i~=j
%             if isSoft
%                 K{i,j}(abs(K{i,j})< 0.01 *maxNorm) = 0;                       
%             else
%                 if A(j,i)==0
%                     K{i,j} = zeros(3);
%                 end
%             end
%         end
% 
%         K_ijMax = max(abs(K{i,j}(:)));
%         K{i,j}(abs(K{i,j})<0.01*K_ijMax) = 0;
% 
%     end
% end

for i=1:1:numOfDGs
    for j=1:1:numOfDGs
        if isSoft
            K{i,j}(abs(K{i,j}) < 0.1 * maxNorm) = 0;
        else
            if A(j,i) == 0
                K{i,j} = zeros(3);
            end
        end

        K_ijMax = max(abs(K{i,j}(:)));
        K{i,j}(abs(K{i,j}) < 0.1 * K_ijMax) = 0;

    end
end



% Load the K_ij values to the DGs
for i = 1:1:numOfDGs
    for j = 1:1:numOfDGs
        DG{i}.K{j} = K{i,j};  
        % K = K{i,j};
    end
end

% Load the K_ij values to the DGs
for i = 1:numOfDGs
    DG{i}.K = cell(1,numOfDGs); % Initialize K as a cell array
    for j = 1:numOfDGs
        DG{i}.K{j} = K(i,j); % Use parentheses for numeric array indexing
    end
end


%% Debugging

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
            W_val = value(constraintMats{idx})
            for j = 1:size(W_val, 1)
                submatrix = W_val(1:j, 1:j);  % Extract principal minor
                if det(submatrix) < 0
                    disp(['Principal minor ', num2str(j), ' is not positive semi-definite.']);
                end
            end
        else
            disp(['Constraint "', constraintTags{idx}, '" is satisfied by ',num2str(feasibility(idx)),' .']);
            W_val = value(constraintMats{idx})
        end
    end

    % costFun0Val = value(costFun0)
    % costFunVal = value(costFun)
    
    Mat_1Eigs = eig(value(Mat_1))
    % Mat_11Val = value(Mat_11);
    % Mat_11Eigs = eig(Mat_11Val)
    % 
    
    Mat_2Eigs = eig(value(Mat_2))
    % Mat_22Val = value(Mat_22);
    % Mat_22Eigs = eig(Mat_22Val)
    % 
    
    Mat_3Eigs = eig(value(Mat_3))
    % Mat_33Val = value(Mat_33);
    % Mat_33Eigs = eig(Mat_33Val)
    % 
    
    Mat_4Eigs = eig(value(Mat_4))
    Mat_44Val = value(Mat_44);
    Mat_44Eigs = eig(Mat_44Val)
    Mat_4_Test1Eigs = eig(value(Mat_4_Test1))
    Mat_4_Test2Eigs = eig(value(Mat_4_Test2))
    Mat_4_Test3Eigs = eig(value(Mat_4_Test3))
    value(P_i)
    value(P_l)
    GammaMat
    value(gammaTilde)
    % 
    
    Mat_5Eigs = eig(value(Mat_5))
    Mat_55Val = value(Mat_55);
    Mat_55Eigs = eig(Mat_55Val)
    Mat_5_Test1Eigs = eig(value(Mat_5_Test1))
    Mat_5_Test2Eigs = eig(value(Mat_5_Test2))
    % 
    
    Mat_6Eigs = eig(value(Mat_6))
    Mat_66Val = value(Mat_66);
    Mat_66Eigs = eig(Mat_66Val)
    Mat_6_Test1Eigs = eig(value(Mat_6_Test1))
    Mat_6_Test2Eigs = eig(value(Mat_6_Test2))


end
end