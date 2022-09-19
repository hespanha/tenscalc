function localVariables_=parameters4minmax(localVariables_)
% Declare input parameters common to the 2 Tenscalc functions:
%   cmex2minmaxCS.m
%   class2minmaxCS.m
%
% This file is part of Tenscalc.
%
% Copyright (C) 2012-22 The Regents of the University of California
% (author: Dr. Joao Hespanha, Raphael Chinchilla).  All rights reserved.

    %% Declare input parameters

    declareParameter(...
        'VariableName','objective',...
        'Description',{
            'Scalar Tcalculus symbolic object to be optimized.'
                      });
    declareParameter(...
        'VariableName','minOptimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be minimized.'
                      });
    declareParameter(...
        'VariableName','minConstraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints for the minimizer. Both equality and inequality'
            'constraints are allowed.'
            ' '
            'ATTENTION: These constrainsts should NOT depend on'
            '           |maxOptimizationVariables|'
            'TO_BE_FIXED: CODE SHOULD CHECK FOR THIS BY COMPUTING DERIVATIVES'
                      });

    declareParameter(...
        'VariableName','maxOptimizationVariables',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'variables to be maximized.'
                      });
    declareParameter(...
        'VariableName','maxConstraints',...
        'DefaultValue',{},...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the'
            'constraints for the maximizer. Both equality and inequality'
            'constraints are acceptable.'
            ' '
            'These constraints may depend on both |minOptimizationVariables|'
            'and |maxOptimizationVariables|.'
                      });

    declareParameter(...
        'VariableName','outputExpressions',...
        'Description',{
            'Cell-array of Tcalculus symbolic objects representing the '
            'variables to be returned.'
            ' ';
            'TO_BE_FIXED'
            ' ';
            'The following Tcalculus symbolic variables are assigned special values';
            'and can be using in outputExpressions';
            '* |P1lambda1_|,|P1lambda2_|,... |P2lambda1_|,|P2lambda2_|,...'
            '      - Lagrangian multipliers associated with the inequalities constraints';
            '        for player 1 and 2 (in the order that they appear and with the same size';
            '        as the corresponding constraints).';
            '* |P1nu1_|,|P1nu2_|,... ,|P2nu1_|,|P2nu2_|,...'
            '* |P1xnu1_|,|P1xnu2_|,... ,|P2xnu1_|,|P2xnu2_|,...'
            '      - Lagrangian multipliers associated with the equality constraints';
            '        for player 1 and 2 (in the order that they appear and with the same size'
            '        as the corresponding constraints). The P1x and P2x variables correspond';
            '        to the latentConstraints.'
            '* |Hess_| - Hessian matrix used by the (last) Newton step to update';
            '            the primal variables (not including |addEye2Hessian|).'
            ' ';
            'ATTENTION: To be able to include these variables as input parameters,';
            '           they have to be previously created outside *with the appropriate sizes*.'
            '           Eventually, their values will be overridden by the solver'
            '           to reflect the values listed above.'
                      });

    declareParameter(...
        'VariableName','addEye2Hessian',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true|, adds to the Newton matrix identity matrices scaled by small constants.';
            ' ';
            'One scaled identity matrix equal to';
            '           addEye2Hessian1 * eye(# primal variables)'
            'is added to the matrix of 2nd derivatives of the Lagrangian (Hessian).'
            ' '
            'A seconds scaled identity matrix equal to';
            '           addEye2Hessian2 * eye(# equality constraints)'
            'is added to the diagonal block of the Newton matrix that corresponds to the equality'
            'constraints, which makes factorization of the Newton matrix numerically more stable.'
            ' '
            'Both effects improve the robustness of the solver, but may lead to slower convergence.';
            ' ';
            'The constants |addEye2Hessian1| and |addEye2Hessian2| can be set at solve time'
            'through input parameters to the solve function.'
            ' '
            'A typical choices for |addEye2Hessian1| and |addEye2Hessian2| is the square root'
            'of the machine precision.'
            ' '
            'The values of |addEye2Hessian1| and |addEye2Hessian2| can be viewed by setting'
            '|solverVerboseLevel| to 3.'
                      });

    declareParameter(...
        'VariableName','adjustAddEye2Hessian',...
        'DefaultValue',false,...
        'AdmissibleValues',{false,true},...
        'Description',{
            'When |true|, the values of the parameters |addEye2Hessian1| and |addEye2Hessian2|';
            'are adjusted in real-time by the solver.';
            ' ';
            'ATTENTION: currently addEye2Hessian1 is not adjusted';
                      });

    declareParameter(...
        'VariableName','addEye2Hessian1tolerance',...
        'DefaultValue',1e-6,...
        'Description',{
            'When |adjustAddEye2Hessian|=|true|, waits until |addEye2Hessian1| becomes smaller than this value.';
            ' ';
            'This parameter is ignored when |adjustAddEye2Hessian|=|false|.'
            ' ';
            'ATTENTION: currently addEye2Hessian1 is not adjusted,'
            '           so addEye2Hessian1 must be initialized below this value.';
                      });

end


% LocalWords:  localVariables tenscalc cmex equilibriumLatentCS Joao
% LocalWords:  declareParameter VariableName Tcalculus DefaultValue
% LocalWords:  optimizationVariables latentVariables xnu addEye nd
% LocalWords:  latentConstraints outputExpressions AdmissibleValues
% LocalWords:  adjustAddEye LDL useLDL solverVerboseLevel
