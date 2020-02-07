function [isRecurrent, affineControllers] = testAffineRecurrencePiecewise(sys,initTargetSet, partition, horizon)
    % testAffineRecurrence determines if there is an affine controller that
    % drives the system from the initial/target set back into itself over
    % the horizon while satisfying all constraints.
    
    assert(horizon <= sys.T);
    assert(horizon > 0);

    % create a copy of the system to add target constraint
    csys = copy(sys);
    csys.Xterm = initTargetSet;
    
    yalmipOptions = sdpsettings('verbose', 1, 'solver', ''); % options for the LP solver
    
    affineControllers = cell(0);
    
    % bisect?
    isRecurrent = 1;
    for i = 1:length(partition(:))
        
        if isEmptySet(partition(i))
            continue;
        end
        
        initialConditions = computeInitialConditions(csys, horizon, partition(i));
        admissibleTrajectories = computeLiftedAdmissibleTrajectories(csys, horizon);

        [diagnostics, affCont] = computeAffineController(csys, initialConditions, admissibleTrajectories, yalmipOptions);
        affineControllers{end+1} = affCont;
        
        if diagnostics.problem ~= 0
            isRecurrent = 0;
            if diagnostics.problem ~= 1
                disp("YALMIP ERROR: ");
                disp(yalmiperror(diagnostics.problem));
            end
            break
        end
    end
    
end

