#pragma once

// Local Includes
#include "GeneralRayKernel.h"

/**
 * A RayKernel that tracks a ray
 */
class TrackingKernel : public GeneralRayKernel
{
public:
  TrackingKernel(const InputParameters & params);

  static InputParameters validParams();

  virtual void onSegment() override;

  virtual void preTrace() override;

  virtual void postTrace() override;

protected:
  Point sampleFluidVelocityField();

private:
  //velocity field of the fluid
  const std::vector<VariableName> _fluid_velocity;

  //The variable number of the fluid velocity in the system
  std::vector<unsigned int> _fluid_velocity_var_num;


  //particle position
  Point _r;

  //particle velocity
  Point _v;

  //time
  Real _t;

  //particle marching time step
  Real _dt;

  //flag to start marching particle
  bool _particle_should_march = true;

  //output vector containing particle position at each time step
  std::vector<std::vector<Real>> _particle_history;
};