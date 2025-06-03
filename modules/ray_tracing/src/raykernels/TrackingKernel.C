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

TrackingKernel::TrackingKernel(const InputParameters & params) : 
GeneralRayKernel(params),
_fluid_velocity(getVectorVar("fluid_velocity", 0)),
_dt(getParam<Real>("dt"))
{}

void
TrackingKernel::preTrace()
{
  _t = 0;
  _r = currentRay()->currentPoint(); 
  _v = sampleFluidVelocityField();  
  _particle_history.push_back({_t, _r(0), _r(1), _r(2)});
  _particle_should_march = true;
}

void
TrackingKernel::onSegment()
{  
  if (_particle_should_march)
  {
    _t += _dt;
    _r += _v * _dt;
    _particle_history.push_back({_t, _r(0), _r(1), _r(2)});
    _particle_should_march = false;
  }

  if (currentRay()->currentElem()->contains_point(_r))
  {
    _v = sampleFluidVelocityField();
    changeRayStartDirection(_r, _v);
    _particle_should_march = true;
  }
}

void
TrackingKernel::postTrace()
{
  std::ofstream output_file("/home/elkamash/projects/moose/modules/ray_tracing/test_cases/particle_position.csv");
  for (const auto & row : _particle_history)
    output_file << row[0] << ',' << row[1] << ',' << row[2] << ',' << row[3] << '\n';
  output_file.close();
}

Point TrackingKernel::sampleFluidVelocityField()
{
  Moose::ElemPointArg elem_pt_arg = {currentRay()->currentElem(), _r, true};
  Moose::StateArg state_arg(static_cast<unsigned int>(0));
  VectorValue v_f = (*_fluid_velocity)(elem_pt_arg, state_arg);
  return Point(MetaPhysicL::raw_value(v_f(0)),  
               MetaPhysicL::raw_value(v_f(1)),  
               MetaPhysicL::raw_value(v_f(2)));
}

