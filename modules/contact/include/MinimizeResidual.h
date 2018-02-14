//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MINIMIZERESIDUAL_H
#define MINIMIZERESIDUAL_H

#include "ConsoleStreamInterface.h"
#include "ContactLineSearch.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/parallel_object.h"
#include <petscsnes.h>

using namespace libMesh;

class FEProblemBase;

class MinimizeResidual : public ContactLineSearch
{
public:
  MinimizeResidual(FEProblemBase & fe_problem,
                   Real cutback_factor,
                   Real growth_factor,
                   MooseApp & app);

  virtual void linesearch(SNESLineSearch linesearch) override;

  /// The current contact state
  virtual std::set<dof_id_type> * contact_state() override { return &_current_contact_state; }
  /// The old contact state
  virtual std::set<dof_id_type> * old_contact_state() override { return &_old_contact_state; }

  void printContactInfo();

protected:
  /// The entire contact set
  std::set<dof_id_type> _current_contact_state;
  std::set<dof_id_type> _old_contact_state;

  PetscReal _user_ksp_rtol;
  bool _user_ksp_rtol_set;
};

#endif
