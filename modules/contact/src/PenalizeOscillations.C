//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PenalizeOscillations.h"
#include "FEProblemBase.h"
#include "libmesh/petsc_solver_exception.h"
#include "petsc/private/linesearchimpl.h"

PenalizeOscillations::PenalizeOscillations(FEProblemBase & fe_problem,
                                           Real cutback_factor,
                                           Real growth_factor,
                                           MooseApp & app)
  : ContactLineSearch(fe_problem, cutback_factor, growth_factor, app)
{
}

void
PenalizeOscillations::linesearch(SNESLineSearch linesearch)
{
  PetscBool changed_y = PETSC_FALSE, changed_w = PETSC_FALSE;
  PetscErrorCode ierr;
  Vec X, F, Y, W, G;
  SNES snes;
  PetscReal fnorm, xnorm, ynorm, gnorm;
  PetscBool domainerror;

  ierr = SNESLineSearchGetVecs(linesearch, &X, &F, &Y, &W, &G);
  LIBMESH_CHKERR(ierr);
  ierr = SNESLineSearchGetNorms(linesearch, &xnorm, &fnorm, &ynorm);
  LIBMESH_CHKERR(ierr);
  ierr = SNESLineSearchGetSNES(linesearch, &snes);
  LIBMESH_CHKERR(ierr);
  ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_SUCCEEDED);
  LIBMESH_CHKERR(ierr);

  /* precheck */
  ierr = SNESLineSearchPreCheck(linesearch, X, Y, &changed_y);
  LIBMESH_CHKERR(ierr);

  /* temporary update */
  _contact_lambda = 1.;
  ierr = VecWAXPY(W, -_contact_lambda, Y, X);
  LIBMESH_CHKERR(ierr);

  /* compute residual to determine whether contact state has changed since the last non-linear
   * residual evaluation */
  _current_contact_state.clear();
  _newly_captured_nodes.clear();
  _newly_released_nodes.clear();
  ierr = (*linesearch->ops->snesfunc)(snes, W, F);
  LIBMESH_CHKERR(ierr);
  ierr = SNESGetFunctionDomainError(snes, &domainerror);
  LIBMESH_CHKERR(ierr);
  if (domainerror)
  {
    ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN);
    LIBMESH_CHKERR(ierr);
  }
  ierr = VecNorm(F, NORM_2, &gnorm);
  LIBMESH_CHKERR(ierr);

  _communicator.set_union(_current_contact_state);
  _console << "\n";
  if (!_current_contact_state.empty())
  {
    _console << "Node ids in contact: ";
    for (auto & node_id : _current_contact_state)
      _console << node_id << " ";
    _console << "\n";
  }
  else
    _console << "No nodes in contact\n";

  std::set_difference(_current_contact_state.begin(),
                      _current_contact_state.end(),
                      _old_contact_state.begin(),
                      _old_contact_state.end(),
                      std::inserter(_newly_captured_nodes, _newly_captured_nodes.begin()));
  std::set_difference(_old_contact_state.begin(),
                      _old_contact_state.end(),
                      _current_contact_state.begin(),
                      _current_contact_state.end(),
                      std::inserter(_newly_released_nodes, _newly_released_nodes.begin()));

  std::set<dof_id_type> in_then_out;
  std::set<dof_id_type> out_then_in;
  std::set_intersection(_previously_captured_nodes.begin(),
                        _previously_captured_nodes.end(),
                        _newly_released_nodes.begin(),
                        _newly_released_nodes.end(),
                        std::inserter(in_then_out, in_then_out.begin()));
  std::set_intersection(_previously_released_nodes.begin(),
                        _previously_released_nodes.end(),
                        _newly_captured_nodes.begin(),
                        _newly_captured_nodes.end(),
                        std::inserter(out_then_in, out_then_in.begin()));

  while ((!in_then_out.empty() || !out_then_in.empty()) && gnorm > fnorm)
  {
    _console << "Lambda = " << _contact_lambda << "\n";
    _console << "Number of oscillating nodes = " << in_then_out.size() + out_then_in.size() << "\n";
    _console << "Oscillating nodes : ";
    for (auto & node_id : in_then_out)
      _console << node_id << " ";
    for (auto & node_id : out_then_in)
      _console << node_id << " ";
    _console << "\n\n";

    _contact_lambda *= 0.5;
    /* update */
    ierr = VecWAXPY(W, -_contact_lambda, Y, X);
    LIBMESH_CHKERR(ierr);

    _current_contact_state.clear();
    _newly_captured_nodes.clear();
    _newly_released_nodes.clear();
    in_then_out.clear();
    out_then_in.clear();
    ierr = (*linesearch->ops->snesfunc)(snes, W, F);
    LIBMESH_CHKERR(ierr);
    ierr = SNESGetFunctionDomainError(snes, &domainerror);
    LIBMESH_CHKERR(ierr);
    if (domainerror)
    {
      ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN);
      LIBMESH_CHKERR(ierr);
    }
    ierr = VecNorm(F, NORM_2, &gnorm);
    LIBMESH_CHKERR(ierr);

    _communicator.set_union(_current_contact_state);
    if (!_current_contact_state.empty())
    {
      _console << "Node ids in contact: ";
      for (auto & node_id : _current_contact_state)
        _console << node_id << " ";
      _console << "\n";
    }

    std::set_difference(_current_contact_state.begin(),
                        _current_contact_state.end(),
                        _old_contact_state.begin(),
                        _old_contact_state.end(),
                        std::inserter(_newly_captured_nodes, _newly_captured_nodes.begin()));
    std::set_difference(_old_contact_state.begin(),
                        _old_contact_state.end(),
                        _current_contact_state.begin(),
                        _current_contact_state.end(),
                        std::inserter(_newly_released_nodes, _newly_released_nodes.begin()));
    std::set_intersection(_previously_captured_nodes.begin(),
                          _previously_captured_nodes.end(),
                          _newly_released_nodes.begin(),
                          _newly_released_nodes.end(),
                          std::inserter(in_then_out, in_then_out.begin()));
    std::set_intersection(_previously_released_nodes.begin(),
                          _previously_released_nodes.end(),
                          _newly_captured_nodes.begin(),
                          _newly_captured_nodes.end(),
                          std::inserter(out_then_in, out_then_in.begin()));
  }
  _console << "Lambda = " << _contact_lambda << "\n";
  if (in_then_out.empty() && out_then_in.empty())
    _console << "No oscillating nodes\n\n";
  else
  {
    _console << "Number of oscillating nodes = " << in_then_out.size() + out_then_in.size() << "\n";
    _console << "Oscillating nodes : ";
    for (auto & node_id : in_then_out)
      _console << node_id << " ";
    for (auto & node_id : out_then_in)
      _console << node_id << " ";
    _console << "\nHowever, norm reduced so moving on.\n\n";
  }

  ierr = VecScale(Y, _contact_lambda);
  LIBMESH_CHKERR(ierr);
  /* postcheck */
  ierr = SNESLineSearchPostCheck(linesearch, X, Y, W, &changed_y, &changed_w);
  LIBMESH_CHKERR(ierr);

  if (changed_y)
  {
    ierr = VecWAXPY(W, -1., Y, X);
    LIBMESH_CHKERR(ierr);
  }

  if (changed_w || changed_y)
  {
    _current_contact_state.clear();
    _newly_captured_nodes.clear();
    _newly_released_nodes.clear();
    ierr = (*linesearch->ops->snesfunc)(snes, W, F);
    LIBMESH_CHKERR(ierr);
    ierr = SNESGetFunctionDomainError(snes, &domainerror);
    LIBMESH_CHKERR(ierr);
    if (domainerror)
    {
      ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN);
      LIBMESH_CHKERR(ierr);
    }
    _communicator.set_union(_current_contact_state);
    if (!_current_contact_state.empty())
    {
      _console << "Node ids in contact: ";
      for (auto & node_id : _current_contact_state)
        _console << node_id << " ";
      _console << "\n";
    }
    std::set_difference(_current_contact_state.begin(),
                        _current_contact_state.end(),
                        _old_contact_state.begin(),
                        _old_contact_state.end(),
                        std::inserter(_newly_captured_nodes, _newly_captured_nodes.begin()));
    std::set_difference(_old_contact_state.begin(),
                        _old_contact_state.end(),
                        _current_contact_state.begin(),
                        _current_contact_state.end(),
                        std::inserter(_newly_released_nodes, _newly_released_nodes.begin()));
  }

  ierr = VecNorm(Y, NORM_2, &linesearch->ynorm);
  LIBMESH_CHKERR(ierr);
  ierr = VecNorm(W, NORM_2, &linesearch->xnorm);
  LIBMESH_CHKERR(ierr);
  ierr = VecNorm(F, NORM_2, &linesearch->fnorm);
  LIBMESH_CHKERR(ierr);

  /* copy the solution over */
  ierr = VecCopy(W, X);
  LIBMESH_CHKERR(ierr);

  _previously_captured_nodes = _newly_captured_nodes;
  _previously_released_nodes = _newly_released_nodes;
  _old_contact_state = _current_contact_state;
}
