#ifndef _QUADRIC_ERROR_METRIC_H_
#define _QUADRIC_ERROR_METRIC_H_

#include "Base/EigenInc.h"
#include "Base/OpenMeshInc.h"

namespace SoftGL {

struct edge_collapse_structure {
	MyMesh::HalfedgeHandle hf;
	MyMesh::Point np;
	MyMesh::VertexHandle vto;
	MyMesh::VertexHandle vfrom;

	int vto_flag = 0;
	int vfrom_flag = 0;
	Eigen::Matrix4f Q_new;
	float cost;

	bool operator<(const edge_collapse_structure& rhs) const {
		return cost > rhs.cost;
	}
};

void simplification(MyMesh& mesh, float simpRatio);
} // namespace SoftGL

#endif // _QUADRIC_ERROR_METRIC_H_