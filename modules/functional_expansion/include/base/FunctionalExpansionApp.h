/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FUNCTIONAL_EXPANSIONAPP_H
#define FUNCTIONAL_EXPANSIONAPP_H

#include "MooseApp.h"

class FunctionalExpansionApp;

template <>
InputParameters validParams<FunctionalExpansionApp>();

class FunctionalExpansionApp : public MooseApp
{
public:
  FunctionalExpansionApp(InputParameters parameters);
  virtual ~FunctionalExpansionApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* FUNCTIONAL_EXPANSIONAPP_H */
