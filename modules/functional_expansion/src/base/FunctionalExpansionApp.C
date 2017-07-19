/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FunctionalExpansionApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<FunctionalExpansionApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

FunctionalExpansionApp::FunctionalExpansionApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  FunctionalExpansionApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  FunctionalExpansionApp::associateSyntax(_syntax, _action_factory);
}

FunctionalExpansionApp::~FunctionalExpansionApp() {}

// External entry point for dynamic application loading
extern "C" void
FunctionalExpansionApp__registerApps()
{
  FunctionalExpansionApp::registerApps();
}
void
FunctionalExpansionApp::registerApps()
{
  registerApp(FunctionalExpansionApp);
}

// External entry point for dynamic object registration
extern "C" void
FunctionalExpansionApp__registerObjects(Factory & factory)
{
  FunctionalExpansionApp::registerObjects(factory);
}
void
FunctionalExpansionApp::registerObjects(Factory & factory)
{
}

// External entry point for dynamic syntax association
extern "C" void
FunctionalExpansionApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  FunctionalExpansionApp::associateSyntax(syntax, action_factory);
}
void
FunctionalExpansionApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
