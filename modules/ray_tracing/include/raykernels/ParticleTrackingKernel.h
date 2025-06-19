#pragma once

// Local Includes
#include "GeneralRayKernel.h"

/**
 * A RayKernel that tracks a ray
 */
class ParticleTrackingKernel : public GeneralRayKernel
{
public:
  ParticleTrackingKernel(const InputParameters & params);

  static InputParameters validParams();

  virtual void onSegment() override;

  virtual void preTrace() override;

  virtual void postTrace() override;

protected:
  void sampleFluidVariables();

  Point buoyancy();

  Point drag();

  Real dragCoefficient(Real Re_r);

private:
  //particle position
  Point _r;

  //particle velocity
  Point _v;

  //particle diameter
  Real _d;

  //particle density
  Real _rho;

  //particle mass
  Real _m;

  //time
  Real _t;

  //particle marching time step
  Real _dt;

  //remaining time step
  Real _particle_dt;

  //flag to start marching particle
  bool _beginning_time_step;

  //whether to correct the velocity at the element boundaries or not
  bool _ray_refraction;

  //first time step
  bool _first_time_step;

  //velocity field of the fluid
  const std::vector<VariableName> _fluid_velocity;

//The variable number of the fluid velocity in the system
  std::vector<unsigned int> _fluid_velocity_var_num;

  //density field of the fluid
  const VariableName _fluid_density;

  //the variable number of the fluid density in the system
  unsigned int _fluid_density_var_num;

  //fluid velocity value at the particle location
  Point _v_f;

  //fluid density value at the particle location
  Real _rho_f;

  //fluid viscosity
  Real _mu;

  //gravity vector
  Point _g;

  //the total force at the particle position
  Point _F;

  //output vector containing particle position at each time step
  std::vector<std::vector<Real>> _particle_history;
};