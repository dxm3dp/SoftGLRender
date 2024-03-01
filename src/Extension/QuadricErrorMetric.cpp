#include "QuadricErrorMetric.h"
#include <queue>

namespace SoftGL {

void simplification(MyMesh& mesh, float simpRatio) {
	assert(simpRatio >= 0 && simpRatio <= 1);

	int it_num = (1.0f - simpRatio) * mesh.n_vertices();

	// Q matrices
	auto Q = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, Eigen::Matrix4f>(mesh);
	// vertex coordinates
	auto v = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, Eigen::Vector4f>(mesh);
	// vertex delete flag
	auto flag = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, int>(mesh);
	// plane equation coefficients
	auto p = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, Eigen::Vector4f>(mesh);

	// iterate over all faces of the mesh to compute the plane equation
	// plane equation : ax + by + cz + d = 0 (a2 + b2 + c2 = 1)
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
		float a, b, c, d;
		a = mesh.normal(*f_it)[0];
		b = mesh.normal(*f_it)[1];
		c = mesh.normal(*f_it)[2];
		MyMesh::Point tp = mesh.point(f_it->halfedge().to());
		d = -(a * tp[0] + b * tp[1] + c * tp[2]);
		p[*f_it][0] = a;
		p[*f_it][1] = b;
		p[*f_it][2] = c;
		p[*f_it][3] = d;
	}

	// iterate over all vertices of the mesh to compute the Q matrix
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		Eigen::Matrix4f mat;
		mat.setZero();
		// iterate over all adjacent faces of the vertex
		for (MyMesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); ++vf_it) {
			mat += p[*vf_it] * p[*vf_it].transpose();
		}
		Q[*v_it] = mat;
		v[*v_it][0] = mesh.point(*v_it)[0];
		v[*v_it][1] = mesh.point(*v_it)[1];
		v[*v_it][2] = mesh.point(*v_it)[2];
		v[*v_it][3] = 1.0f;
		flag[*v_it] = 0;
	}
	// select all valid pairs(now only vertices in an edge are considered)
	std::priority_queue<edge_collapse_structure, std::vector<edge_collapse_structure>, std::less<edge_collapse_structure>> heap;
	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
		Eigen::Matrix4f newQ = Q[e_it->v0()] + Q[e_it->v1()];
		Eigen::Matrix4f tQ = newQ;
		tQ(3, 0) = 0.0f;
		tQ(3, 1) = 0.0f;
		tQ(3, 2) = 0.0f;
		tQ(3, 3) = 1.0f;
		// compute LU decomposition with full pivoting for deciding if the matrix is invertible
		Eigen::FullPivLU<Eigen::Matrix4f> lu(tQ);
		Eigen::Vector4f newV;
		Eigen::Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
		if (lu.isInvertible()) {
			newV = tQ.inverse() * b;
		} else {
			newV = (v[e_it->v0()] + v[e_it->v1()]) / 2.0f;
		}
		// std::cout << newV << std::endl;
		edge_collapse_structure ts;
		ts.hf = e_it->halfedge(0);
		ts.cost = newV.transpose() * newQ * newV;
		ts.np = {newV[0], newV[1], newV[2]};
		ts.vto = e_it->halfedge(0).to();
		ts.vfrom = e_it->halfedge(0).from();
		ts.Q_new = newQ;
		heap.push(ts);
	}
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	int i = 0;
	while (i < it_num) {
		edge_collapse_structure s = heap.top();
		heap.pop();
		if (mesh.status(s.vfrom).deleted() || mesh.status(s.vto).deleted())
			continue;
		if (s.vto_flag != flag[s.vto] || s.vfrom_flag != flag[s.vfrom])
			continue;

		MyMesh::VertexHandle tvh;
		if (mesh.is_collapse_ok(s.hf)) {
			mesh.collapse(s.hf);
			tvh = s.vto;
			++flag[s.vto];
			++flag[s.vfrom];
		} else if (mesh.is_collapse_ok(mesh.opposite_halfedge_handle(s.hf))) {
			mesh.collapse(mesh.opposite_halfedge_handle(s.hf));
			tvh = s.vfrom;
			++flag[s.vto];
			++flag[s.vfrom];
		} else {
			continue;
		}
		// update Q matrix and vertex coordinates
		mesh.set_point(tvh, s.np);
		Q[tvh] = s.Q_new;
		v[tvh][0] = s.np[0];
		v[tvh][1] = s.np[1];
		v[tvh][2] = s.np[2];
		v[tvh][3] = 1.0f;
		// iterate over all outgoing half edge to update each pair info
		for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(tvh); voh_it.is_valid(); ++voh_it) {
			MyMesh::VertexHandle tt = voh_it->to();
			Eigen::Matrix4f newQ = s.Q_new + Q[tt];
			Eigen::Matrix4f tQ = newQ;
			tQ(3, 0) = 0.0f;
			tQ(3, 1) = 0.0f;
			tQ(3, 2) = 0.0f;
			tQ(3, 3) = 1.0f;
			Eigen::FullPivLU<Eigen::Matrix4f> lu(tQ);
			Eigen::Vector4f newV;
			Eigen::Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
			if (lu.isInvertible()) {
				newV = tQ.inverse() * b;
			} else {
				newV = (v[tvh] + v[tt]) / 2.0f;
			}
			// std::cout << newV << std::endl;
			edge_collapse_structure ts;
			ts.hf = *voh_it;
			ts.cost = newV.transpose() * newQ * newV;
			MyMesh::Point np(newV[0], newV[1], newV[2]);
			ts.np = np;
			ts.vto = tt;
			ts.vto_flag = flag[tt];
			ts.vfrom = tvh;
			ts.vfrom_flag = flag[tvh];
			ts.Q_new = newQ;
			heap.push(ts);
		}
		++i;
	}
	mesh.garbage_collection();
}
} // namespace SoftGL