#pragma once

#include "Base/GLMInc.h"
#include "Viewer/Model.h"

#include <memory>

namespace SoftGL {
struct PairContraction {
  float cost;
  int i, j;
  glm::vec3 v_bar;
  bool operator<(const PairContraction &o) const { return cost > o.cost; }
};

class QuadricErrorMetrics {
public:
  void init(View::Model &model);
};

} // namespace SoftGL