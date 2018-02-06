//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PenalizeAnyChange.h"
#include "FEProblemBase.h"
#include "libmesh/petsc_solver_exception.h"
#include "petsc/private/linesearchimpl.h"

PenalizeAnyChange::PenalizeAnyChange(FEProblemBase & fe_problem,
                                     Real cutback_factor,
                                     Real growth_factor,
                                     MooseApp & app)
  : ContactLineSearch(fe_problem, cutback_factor, growth_factor, app)
{
}

void
PenalizeAnyChange::linesearch(SNESLineSearch linesearch)
{
  PetscBool changed_y = PETSC_FALSE, changed_w = PETSC_FALSE;
  bool need_contact_eval = false;
  PetscErrorCode ierr;
  Vec X, F, Y, W, G;
  SNES snes;
  PetscReal gnorm, xnorm, ynorm;
  PetscBool domainerror;

  ierr = SNESLineSearchGetVecs(linesearch, &X, &F, &Y, &W, &G);
  LIBMESH_CHKERR(ierr);
  ierr = SNESLineSearchGetNorms(linesearch, &xnorm, &gnorm, &ynorm);
  LIBMESH_CHKERR(ierr);
  ierr = SNESLineSearchGetSNES(linesearch, &snes);
  LIBMESH_CHKERR(ierr);
  ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_SUCCEEDED);
  LIBMESH_CHKERR(ierr);

  /* precheck */
  ierr = SNESLineSearchPreCheck(linesearch, X, Y, &changed_y);
  LIBMESH_CHKERR(ierr);

  /* temporary update */
  ierr = VecWAXPY(W, -_contact_lambda, Y, X);
  LIBMESH_CHKERR(ierr);

  /* compute residual to determine whether contact state has changed since the last non-linear
   * residual evaluation */
  _current_contact_state.clear();
  ierr = (*linesearch->ops->snesfunc)(snes, W, F);
  LIBMESH_CHKERR(ierr);
  ierr = SNESGetFunctionDomainError(snes, &domainerror);
  LIBMESH_CHKERR(ierr);
  if (domainerror)
  {
    ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN);
    LIBMESH_CHKERR(ierr);
  }

  _communicator.set_union(_current_contact_state, 0);
  _console << "With initial lambda of " << _contact_lambda << ", " << _current_contact_state.size()
           << " nodes in contact.\n";

  bool contact_changing_this_nl(_current_contact_state != _old_contact_state);
  if (contact_changing_this_nl)
    _contact_changing_this_timestep = true;

  if (_contact_changing_this_timestep)
  {
    _old_contact_lambda = _contact_lambda;
    if (contact_changing_this_nl)
    {
      _contact_lambda = std::max(1e-3, _contact_lambda * _cutback_factor);
      if (fabs(_old_contact_lambda - _contact_lambda) > std::numeric_limits<Real>::epsilon())
      {
        _console << "Changed contact state!!! Decreasing lambda to " << _contact_lambda << ".\n";

        /* update */
        ierr = VecWAXPY(W, -_contact_lambda, Y, X);
        LIBMESH_CHKERR(ierr);
        need_contact_eval = true;
      }
      else
        _console << "Changed contact state!!! But lambda already at lambda_min.\n";
    }
    else
    {
      _contact_lambda = std::min(1., _contact_lambda * _growth_factor);
      if (fabs(_old_contact_lambda - _contact_lambda) > std::numeric_limits<Real>::epsilon())
      {
        _console << "Unchanged contact state. Increasing contact lambda to " << _contact_lambda
                 << ".\n";

        /* update */
        ierr = VecWAXPY(W, -_contact_lambda, Y, X);
        LIBMESH_CHKERR(ierr);

        _current_contact_state.clear();
        ierr = (*linesearch->ops->snesfunc)(snes, W, G);
        LIBMESH_CHKERR(ierr);
        ierr = SNESGetFunctionDomainError(snes, &domainerror);
        LIBMESH_CHKERR(ierr);
        if (domainerror)
        {
          ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN);
          LIBMESH_CHKERR(ierr);
        }
        _communicator.set_union(_current_contact_state, 0);
        if (_current_contact_state == _old_contact_state)
        {
          ierr = VecCopy(G, F);
          LIBMESH_CHKERR(ierr);
        }
        else
        {
          _console << "Increasing contact lambda changed the contact state, so we revert to the "
                      "old value.\n";
          _contact_lambda = _old_contact_lambda;
          _current_contact_state = _old_contact_state;
          // revert working vector
          ierr = VecWAXPY(W, -_contact_lambda, Y, X);
          LIBMESH_CHKERR(ierr);
        }
      }
      else
        _console << "Unchanged contact state. Lambda is already unity.\n";
    }
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

  if (need_contact_eval || changed_w || changed_y)
  {
    _current_contact_state.clear();
    ierr = (*linesearch->ops->snesfunc)(snes, W, F);
    LIBMESH_CHKERR(ierr);
    ierr = SNESGetFunctionDomainError(snes, &domainerror);
    LIBMESH_CHKERR(ierr);
    if (domainerror)
    {
      ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN);
      LIBMESH_CHKERR(ierr);
    }
    _communicator.set_union(_current_contact_state, 0);
    _console << "After modifying lambda, " << _current_contact_state.size()
             << " nodes in contact.\n";
  }

  _old_contact_state = _current_contact_state;

  ierr = VecNormBegin(Y, NORM_2, &linesearch->ynorm);
  LIBMESH_CHKERR(ierr);
  ierr = VecNormBegin(W, NORM_2, &linesearch->xnorm);
  LIBMESH_CHKERR(ierr);
  ierr = VecNormEnd(Y, NORM_2, &linesearch->ynorm);
  LIBMESH_CHKERR(ierr);
  ierr = VecNormEnd(W, NORM_2, &linesearch->xnorm);
  LIBMESH_CHKERR(ierr);

  ierr = VecNorm(F, NORM_2, &linesearch->fnorm);
  LIBMESH_CHKERR(ierr);

  /* copy the solution over */
  ierr = VecCopy(W, X);
  LIBMESH_CHKERR(ierr);
}
