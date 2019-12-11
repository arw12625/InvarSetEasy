function [isRecurrent, affineController] = testAffineRecurrence(sys,initTargetSet,horizon)
    % testAffineRecurrence determines if there is an affine controller that
    % drives the system from the initial/target set back into itself over
    % the horizon while satisfying all constraints.
    
    assert(horizon <= sys.T);
    assert(horizon > 0);

    % create a copy of the system to add target constraint
    csys = copy(sys);
    csys.Xterm = csys.Xterm & initTargetSet;
    
    initialConditions = computeInitialConditions(csys, initTargetSet, horizon);
    admissibleTrajectories = computeLiftedAdmissibleTrajectories(csys, horizon);
    
    yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver
    [diagnostics, affineController] = computeAffineController(csys, initialConditions, admissibleTrajectories, yalmipOptions);
    
    isRecurrent = 1 - diagnostics.problem;
    
end

