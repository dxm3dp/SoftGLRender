#include "QuadricErrorMetric.h"

namespace SoftGL {
void simplification(MyMesh& mesh, float simpRatio) {
	assert(simpRatio >= 0 && simpRatio <= 1);

	int it_num = (1.0f - simpRatio) * mesh.n_vertices();

	// Q
	auto Q = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, float>(mesh);
}
} // namespace SoftGL