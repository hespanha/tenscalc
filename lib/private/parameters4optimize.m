function localVariables_=parameters4optimize(localVariables_)
% Declare input parameters common to the 4 tenscalc functions:
%   cmex2optimizeCS.m
%   class2optimizeCS.m

    %% Declare input parameters
    declareParameter(...
        'VariableName','objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized.'
                      });
    declareParameter(...
        'VariableName','optimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be optimized.'
                      });
    declareParameter(...
        'VariableName','sensitivityVariables',...
        'DefaultValue',{},...        
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'optimization variables with respect to which we want to compute cost';
            'sensitivity.'
                      });
    declareParameter(...
        'VariableName','constraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints. Both equality and inequality constraints'
            'are allowed.'
                      });
    declareParameter(...
        'VariableName','outputExpressions',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'variables to be returned.';
            ' ';
            'The following Tcalculus symbolic variables are assigned special values';
            'and can be using in outputExpressions';
            '* |lambda1_|,|lambda2_|,... - Lagrangian multipliers associated with';
            '                         the inequalities constraints';
            '                         (in the order that they appear and with';
            '                          the same size as the corresponding constraints)';
            '* |nu1_|,|nu2_|,...    - Lagrangian multipliers associated with';
            '                         the equality constraints';
            '                         (in the order that they appear and with';
            '                          the same size as the corresponding constraints)';
            '* |Hess_| - Hessian matrix used by the (last) Newton step.'
            '* |dHess_| - D factor in the LDL factorization of the Hessian matrix'
            '             used by the (last) Newton step.'
            '* |Grad_| - gradient of Lagrangian at the (last) Newton step.'
            '* |mu_|   - barrier parameter at the (last) Newton step.'
            '* |u_|    - vector stacked with all primal variables at the (last) Newton step.'
            '* |F_|    - vector stacked with all equalities at the (last) Newton step.'
            '* |G_|    - vector stacked with all inequalities at the (last) Newton step.'
            '* |nu_|   - vector stacked with all dual equality variables at the (last) Newton step.'
            '* |lambda_|   - vector stacked with all dual inequality variables at the (last) Newton step.'
            ' ';
            'ATTENTION: To be able to include these variables as input parameters,';
            '           they have to be previously created outside *with the appropriate sizes*.'
            '           Eventually, their values will be overridden by the solver'
            '           to reflect the values listted above.'
                      });

    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...        
        'Description',{
            'When |true|, adds to the Newton matrix identity matrices scaled by a small constant.';
            ' ';
            'One scaled identity matrix equal to';
            '           addEye2Hessian1 * eye(# primal variables)'
            'is added to the matrix of 2nd derivatives of the Lagragian (Hessian), helping this'
            'matrix to become positive definite and moving the Newton''s';
            'search direction towards the gradient descent of the Lagragian (and away from the'
            'pure Newton direction).';
            ' ';
            'A seconds scaled identity matrix equal to';
            '           addEye2Hessian2 * eye(# equality constraints)'
            'is added to the diagonal block of the Newton matrix that corresponds to the equality'
            'constraints, turning it slighly negative definite, which makes factrorization'
            'of the Newton matrix numerically more stable.'
            ' '
            'Both effects improve the robustness of the solver, but may lead to slower convergence.';
            ' ';
            'The constants ''addEye2Hessian1'' and ''addEye2Hessian2'' can be set at solve time'
            'through input parameters to the solve function and can also be adjusted by the solver'
            'at run time. See ''adjustAddEye2Hessian''.'
            ' '
            'A typical choices for ''addEye2Hessian1'' and ''addEye2Hessian2'' is the square root'
            'of the machine precision.'
            ' '
            'For non-convex problems, one can try to increase this parameter when';
            'the Newton direction actually causes an increase of the Lagragian.'
                      });
    
    declareParameter(...
        'VariableName','adjustAddEye2Hessian',...
        'DefaultValue',true,...
        'AdmissibleValues',{false,true},...        
        'Description',{
            'When |true|, the values of the parametes ''addEye2Hessian1'' and ''addEye2Hessian2''';
            'are adjusted in rwal-time by the solver.';
            ' '
            'This is only possivle when using LDL factorization (''useLDL'' set to true) and will result';
            'slightly slower solvers. Ideally, one would try a few test runs with ''adjustAddEye2Hessian''';
            'set to ''true'' to learn good values for ''addEye2Hessian1'' and ''addEye2Hessian2'' and then turn';
            '''adjustAddEye2Hessian'' to ''false''.';
            ' ';
            'The values of ''addEye2Hessian1'' and ''addEye2Hessian2'' can be viewed by setting'
            '''solverVerboseLevel'' to 3.'
                      });
    
end

