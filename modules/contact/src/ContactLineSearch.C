//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ContactLineSearch.h"
#include "FEProblemBase.h"
#include "libmesh/petsc_solver_exception.h"
#include "petsc/private/linesearchimpl.h"

ContactLineSearch::ContactLineSearch(FEProblemBase & fe_problem,
                                     Real cutback_factor,
                                     Real growth_factor,
                                     MooseApp & app)
  : ConsoleStreamInterface(app),
    ParallelObject(app),
    _fe_problem(fe_problem),
    _contact_changing_this_timestep(false),
    _cutback_factor(cutback_factor),
    _growth_factor(growth_factor),
    _contact_lambda(1.),
    _old_contact_lambda(1.),
    _nl_its(0)
{
}
