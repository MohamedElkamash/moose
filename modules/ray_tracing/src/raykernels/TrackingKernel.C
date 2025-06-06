#include "TrackingKernel.h"

registerMooseObject("RayTracingApp", TrackingKernel);

InputParameters
TrackingKernel::validParams()
{
  auto params = GeneralRayKernel::validParams();
  params.addClassDescription("A RayKernel that tracks a particle");
  params.addRequiredParam<Real>("dt", "track integration time step");
  params.addRequiredParam<std::vector<VariableName>>("fluid_velocity", "The velocity field of the fluid.");
  return params;
}

TrackingKernel::TrackingKernel(const InputParameters & params) : 
GeneralRayKernel(params),
_fluid_velocity(getParam<std::vector<VariableName>>("fluid_velocity")),
_dt(getParam<Real>("dt"))
{ 
  for (int i=0; i<3; ++i)
    _fluid_velocity_var_num.push_back( 
      _fe_problem.getVariable(_tid, _fluid_velocity[i], Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD).number());
}

void
TrackingKernel::preTrace()
{
  _t = 0;
  _r = currentRay()->currentPoint(); 
  _v = sampleFluidVelocityField();  
  _particle_history.push_back({_t, _r(0), _r(1), _r(2)});
  _beginning_time_step = true;
  _v.print();
}

void
TrackingKernel::onSegment()
{
  if (_beginning_time_step)
  {
    _particle_dt = _dt;
    _beginning_time_step = false;
  }

  Real segment_dt = _current_segment_length / _v.norm();
  
  if (_particle_dt < segment_dt)
  {
    _t += _particle_dt;
    _r += _particle_dt * _v;
    _v = sampleFluidVelocityField();
    changeRayStartDirection(_r, _v);
    _beginning_time_step = true;
    _particle_history.push_back({_t, _r(0), _r(1), _r(2)});
  }
  else
  {
    _t += segment_dt;
    _r = currentRay()->currentPoint();
    _particle_dt -= segment_dt;
  }
}

void
TrackingKernel::postTrace()
{
  //_r = currentRay()->currentPoint();
  //_particle_history.push_back({_t, _r(0), _r(1), _r(2)});
  //std::ofstream output_file("/home/elkamash/projects/moose/modules/ray_tracing/test_cases/particle_position.csv");
    std::ofstream output_file("/Users/elkamm/projects/moose/modules/ray_tracing/test_cases/particle_position.csv");
  for (const auto & row : _particle_history)
    output_file << row[0] << ',' << row[1] << ',' << row[2] << ',' << row[3] << '\n';
  output_file.close();
}

Point TrackingKernel::sampleFluidVelocityField()
{
  std::vector<Real> v_f(3);
  for (int i=0; i<3; ++i)
    v_f[i] = _fe_problem.getSystem(_fluid_velocity[i]).point_value(_fluid_velocity_var_num[i], _r, currentRay()->currentElem());
  return Point(v_f[0], v_f[1], v_f[2]);
}

