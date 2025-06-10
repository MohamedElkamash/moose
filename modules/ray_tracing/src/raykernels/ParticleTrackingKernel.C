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
  params.addRequiredParam<Real>("fluid_viscosity", "The viscosity of the fluid");
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
_mu(getParam<Real>("fluid_viscosity")),
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
    _F = drag();
    _v += _dt * _F / _m;
    _v.print();
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

Point ParticleTrackingKernel::drag()
{
  Real Re_r = _rho_f * _d / _mu * ((_v_f - _v).norm());
  Real C_drag = dragCoefficient(Re_r);
  Real tau = 4*_rho*pow(_d,2) / (3*_mu*C_drag*Re_r);
  return _m * (_v_f - _v) / tau;
}

Real ParticleTrackingKernel::dragCoefficient(Real Re_r)
{
  std::cout << "Re = " << Re_r << '\n';
  if (Re_r <= 0)
    return 0;
  else if (Re_r < 0.01)
    return 24.0/Re_r * (1 + 3.0/16.0*Re_r);
  else 
    return 24.0/Re_r * (1 + 0.1315 * pow(Re_r, 0.82 - 0.05 * log10(Re_r))); 
}