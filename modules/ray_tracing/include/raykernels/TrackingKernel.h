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
  const VectorMooseVariable * _fluid_velocity;

  Point _r;

  Point _v;

  Real _dt;
};