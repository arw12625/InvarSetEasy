function [isRecurrent, affineController] = testAffineRecurrence(sys,initTargetSet,horizon)
    % testAffineRecurrence determines if there is an affine controller that
    % drives the system from the initial/target set back into itself over
    % the horizon while satisfying all constraints.
    
    assert(horizon <= sys.T);
    assert(horizon > 0);

    % create a copy of the system to add target constraint
    csys = copy(sys);
    csys.Xterm = csys.Xterm & initTargetSet;
    
    initialConditions = computeInitialConditions(csys, horizon, initTargetSet);
    admissibleTrajectories = computeLiftedAdmissibleTrajectories(csys, horizon);
    2
    yalmipOptions = sdpsettings('verbose', 2, 'solver', 'gurobi'); % options for the LP solver
    tic
    [diagnostics, affineController] = computeAffineController(csys, initialConditions, admissibleTrajectories, yalmipOptions);
    toc
    diagnostics
    if diagnostics.problem == 0
        isRecurrent = 1;
    else
        isRecurrent = 0;
        if diagnostics.problem ~= 1
            disp("YALMIP ERROR: ");
            disp(yalmiperror(diagnostics.problem));
        end
    end
    
end

