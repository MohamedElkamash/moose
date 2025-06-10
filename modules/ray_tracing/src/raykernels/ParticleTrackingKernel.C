#include "ParticleTrackingKernel.h"

registerMooseObject("RayTracingApp", ParticleTrackingKernel);

InputParameters
ParticleTrackingKernel::validParams()
{
  auto params = GeneralRayKernel::validParams();
  params.addClassDescription("A RayKernel that tracks a particle");
  params.addRequiredParam<Real>("particle_diameter", "diameter of the particle");
  params.addRequiredParam<Real>("particle_density", "density of the particle");
  params.addRequiredParam<Point>("initial_position", "initial position of the particle");
  params.addRequiredParam<Point>("initial_velocity", "initial velocity of the particle");
  params.addRequiredParam<Real>("dt", "track integration time step");
  params.addRequiredParam<std::vector<VariableName>>("fluid_velocity", "The velocity field of the fluid.");
  params.addRequiredParam<VariableName>("fluid_density", "The density field of the fluid.");
  params.addRequiredParam<Point>("gravity", "The gravity vector");
  params.addParam<bool>("ray_refraction", false, "if true, corrects the direction of the particle at each intersection");
  return params;
}

ParticleTrackingKernel::ParticleTrackingKernel(const InputParameters & params) : 
GeneralRayKernel(params),
_r(getParam<Point>("initial_position")),
_v(getParam<Point>("initial_velocity")),
_d(getParam<Real>("particle_diameter")),
_rho(getParam<Real>("particle_density")),
_m(_rho * 3.14159 * _d * _d * _d / 6.0),
_t(0),
_dt(getParam<Real>("dt")),
_ray_refraction(getParam<bool>("ray_refraction")),
_fluid_velocity(getParam<std::vector<VariableName>>("fluid_velocity")),
_fluid_density(getParam<VariableName>("fluid_density")),
_g(getParam<Point>("gravity"))
{ 
  for (int i=0; i<3; ++i)
    _fluid_velocity_var_num.push_back( 
      _fe_problem.getVariable(_tid, _fluid_velocity[i], Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD).number());

  _fluid_density_var_num = 
      _fe_problem.getVariable(_tid, _fluid_density, Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD).number();
}

void
ParticleTrackingKernel::preTrace()
{
  _beginning_time_step = true; 
  _particle_history.push_back({_t, _r(0), _r(1), _r(2)});
}

void
ParticleTrackingKernel::onSegment()
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
    sampleFluidVariables();
    _F = buoyancy();
    _v += _dt * _F / _m;
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
ParticleTrackingKernel::postTrace()
{
  //_r = currentRay()->currentPoint();
  //_particle_history.push_back({_t, _r(0), _r(1), _r(2)});
  //std::ofstream output_file("/home/elkamash/projects/moose/modules/ray_tracing/test_cases/particle_position.csv");
    std::ofstream output_file("/Users/elkamm/projects/moose/modules/ray_tracing/test_cases/particle_position.csv");
  for (const auto & row : _particle_history)
    output_file << row[0] << ',' << row[1] << ',' << row[2] << ',' << row[3] << '\n';
  output_file.close();
}

void ParticleTrackingKernel::sampleFluidVariables()
{
  std::vector<Real> v_f(3);
  for (int i=0; i<3; ++i)
    v_f[i] = _fe_problem.getSystem(_fluid_velocity[i]).point_value(_fluid_velocity_var_num[i], _r, currentRay()->currentElem());
  _v_f = Point(v_f[0], v_f[1], v_f[2]);

  _rho_f = _fe_problem.getSystem(_fluid_density).point_value(_fluid_density_var_num, _r, currentRay()->currentElem());
}

Point ParticleTrackingKernel::buoyancy()
{
  return _m * (1.0 - _rho_f/_rho) * _g;
}

