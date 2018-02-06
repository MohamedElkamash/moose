//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PENALIZEANYCHANGE_H
#define PENALIZEANYCHANGE_H

#include "ConsoleStreamInterface.h"
#include "ContactLineSearch.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/parallel_object.h"
#include <petscsnes.h>

using namespace libMesh;

class FEProblemBase;

class PenalizeAnyChange : public ContactLineSearch
{
public:
  PenalizeAnyChange(FEProblemBase & fe_problem,
                    Real cutback_factor,
                    Real growth_factor,
                    MooseApp & app);

  virtual void linesearch(SNESLineSearch linesearch) override;

  /// The current contact state
  virtual std::set<dof_id_type> * contact_state() override { return &_current_contact_state; }

protected:
  /// The entire contact set
  std::set<dof_id_type> _current_contact_state;
  std::set<dof_id_type> _old_contact_state;
};

#endif
