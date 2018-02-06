//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CONTACTLINESEARCH_H
#define CONTACTLINESEARCH_H

#include "ConsoleStreamInterface.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/parallel_object.h"
#include <petscsnes.h>

using namespace libMesh;

class FEProblemBase;

class ContactLineSearch : public PetscNonlinearSolver<Real>::ComputeLineSearchObject,
                          public ConsoleStreamInterface,
                          public ParallelObject
{
public:
  ContactLineSearch(FEProblemBase & fe_problem,
                    Real cutback_factor,
                    Real growth_factor,
                    MooseApp & app);

  /// Returns a writeable reference to the _contact_changing_this_timestep flag
  bool & contactChangingThisTimestep() { return _contact_changing_this_timestep; }

  /// Our lambda
  Real & lambda() { return _contact_lambda; }

  /// The current contact state
  virtual std::set<dof_id_type> * contact_state() { return nullptr; }

  /// newly captured nodes
  virtual std::set<dof_id_type> * newly_captured_nodes() { return nullptr; }

  /// newly released nodes
  virtual std::set<dof_id_type> * newly_released_nodes() { return nullptr; }

protected:
  FEProblemBase & _fe_problem;

  /// Flag that tracks the timestep status for contact problems. Necessary for custom line search
  bool _contact_changing_this_timestep;

  /// cutback factor
  Real _cutback_factor;

  /// growth factor
  Real _growth_factor;

  /// line search scaling factor for contact
  Real _contact_lambda;

  /// old contact lambda
  Real _old_contact_lambda;
};

#endif
