#include "TrackingKernel.h"

registerMooseObject("RayTracingApp", TrackingKernel);

InputParameters
TrackingKernel::validParams()
{
  auto params = GeneralRayKernel::validParams();
  params.addClassDescription("A RayKernel that tracks a particle");
  params.addRequiredParam<Real>("dt", "track integration time step");
  params.addRequiredCoupledVar("fluid_velocity", "velocity field of the fluid");
  return params;
}

TrackingKernel::TrackingKernel(const InputParameters & params) : GeneralRayKernel(params),
_fluid_velocity(getVectorVar("fluid_velocity", 0)),
_dt(getParam<Real>("dt"))
{}

void
TrackingKernel::preTrace()
{
    _r = currentRay()->currentPoint(); //particle initial position
    Moose::ElemPointArg elem_pt_arg = {currentRay()->currentElem(), _r, true};
    Moose::StateArg state_arg(static_cast<unsigned int>(0));
    VectorValue v_f = (*_fluid_velocity)(elem_pt_arg, state_arg);
    _v = Point(MetaPhysicL::raw_value(v_f(0)),  
               MetaPhysicL::raw_value(v_f(1)),  
               MetaPhysicL::raw_value(v_f(2))); //particle initial velocity
    _v.print();
}

void
TrackingKernel::onSegment()
{
    //_current_segment_start.print();
    //std::cout << "\n";
    // const std::shared_ptr<Ray> ray = currentRay();
    // std::cout << ray->getInfo() << '\n';
}

void
TrackingKernel::postTrace()
{
    //std::cout << currentRay()->getInfo() << '\n';
}

