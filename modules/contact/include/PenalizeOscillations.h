//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PENALIZEOSCILLATIONS_H
#define PENALIZEOSCILLATIONS_H

#include "ConsoleStreamInterface.h"
#include "ContactLineSearch.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/parallel_object.h"
#include <petscsnes.h>

using namespace libMesh;

class FEProblemBase;

class PenalizeOscillations : public ContactLineSearch
{
public:
  PenalizeOscillations(FEProblemBase & fe_problem,
                       Real cutback_factor,
                       Real growth_factor,
                       MooseApp & app);

  virtual void linesearch(SNESLineSearch linesearch) override;

  /// newly captured nodes
  virtual std::set<dof_id_type> * newly_captured_nodes() override { return &_newly_captured_nodes; }

  /// newly released nodes
  virtual std::set<dof_id_type> * newly_released_nodes() override { return &_newly_released_nodes; }

protected:
  /// newly captured nodes
  std::set<dof_id_type> _newly_captured_nodes;
  /// newly released nodes
  std::set<dof_id_type> _newly_released_nodes;

  /// nodes captured on previous nl
  std::set<dof_id_type> _previously_captured_nodes;
  /// nodes released on previous nl
  std::set<dof_id_type> _previously_released_nodes;
};

#endif
