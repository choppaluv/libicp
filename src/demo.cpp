/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libicp.
Authors: Andreas Geiger

libicp is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or any later version.

libicp is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libicp; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA
*/

// Demo program showing how libicp can be used

#include <iostream>
#include <random>
#include <cmath>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/math/matrix44.h>
#include "icpPointToPlane.h"
#include "MeshModel.h"
#include "config.h"

using namespace std;
void visualizeResult(MyMeshOcf& m, kdtree::KDTreeArray pointset, bool pointToPoint = false);
int main(int argc, char** argv) {

	// define a 3 dim problem with 10000 model points
	// and 10000 template points:
	// This configuration can be changed in config.h file.
	int32_t dim = DIM;
	int32_t num = NUM_SAMPLE_POINTS;

	// allocate model and template memory
	double* M = (double*)calloc(3 * num, sizeof(double));
	double* T = (double*)calloc(3 * num, sizeof(double));

	MyMeshOcf m1, m2;
	const char *filename1 = "../source_stl/translated_new.stl";
	const char *filename2 = "../source_stl/translated2_new.stl";

	// import STL file
	int dummy = 0;
	if (vcg::tri::io::ImporterSTL<MyMeshOcf>::Open(m1, filename1, dummy, 0) != 0) {
		std::printf("Error reading file  %s\n", filename1);
		exit(0);
	}

	if (vcg::tri::io::ImporterSTL<MyMeshOcf>::Open(m2, filename2, dummy, 0) != 0) {
		std::printf("Error reading file  %s\n", filename2);
		exit(0);
	}

	// Store all vertices of each mesh.
	kdtree::KDTreeArray A, B;

	for (MyMeshOcf::VertexIterator vi = m1.vert.begin();vi != m1.vert.end();++vi) {
		std::vector<float> point;
		point.push_back((*vi).P()[0]);
		point.push_back((*vi).P()[1]);
		point.push_back((*vi).P()[2]);
		A.push_back(point);
	}

	for (MyMeshOcf::VertexIterator vi = m2.vert.begin();vi != m2.vert.end();++vi) {
		std::vector<float> point;
		point.push_back((*vi).P()[0]);
		point.push_back((*vi).P()[1]);
		point.push_back((*vi).P()[2]);
		B.push_back(point);
	}

	std::random_device randDev;
	std::mt19937 generator(randDev());
	std::uniform_int_distribution<int> distrA(0, m1.VN() - 1);
	std::uniform_int_distribution<int> distrB(0, m2.VN() - 1);

	cout << "dicom size: " << m1.VN() << endl;
	cout << "stl size: " << m2.VN() << endl;
	cout << "sample size: " << num << endl;
	set<int> indices1;
	set<int> indices2;

	// set model and template points
	for (int i = 0; i < num; i++) {
		int idx1 = distrA(generator);
		int idx2 = distrB(generator);
		indices1.insert(idx1);
		indices2.insert(idx2);

		M[i * 3 + 0] = A.at(idx1).at(0);
		M[i * 3 + 1] = A.at(idx1).at(1);
		M[i * 3 + 2] = A.at(idx1).at(2);
		T[i * 3 + 0] = B.at(idx2).at(0);
		T[i * 3 + 1] = B.at(idx2).at(1);
		T[i * 3 + 2] = B.at(idx2).at(2);
	}

	// Color the sampled points
	int i = 0;
	for (MyMeshOcf::VertexIterator vi = m1.vert.begin();vi != m1.vert.end();++vi) {
		if (indices1.find(i) != indices1.end())
			(*vi).SetS();
		else
			(*vi).ClearS();
		i++;
	}

	i = 0;
	for (MyMeshOcf::VertexIterator vi = m2.vert.begin();vi != m2.vert.end();++vi) {
		if (indices2.find(i) != indices2.end())
			(*vi).SetS();
		else
			(*vi).ClearS();
		i++;
	}

	// Since Vertex Coloring is not supported for exporting stl, set the color for the faces.
	for (MyMeshOcf::FaceIterator fi = m1.face.begin(); fi != m1.face.end(); ++fi) if (!(*fi).IsD()) {
		if ((*fi).V(0)->IsS() || (*fi).V(1)->IsS() || (*fi).V(2)->IsS()) {
			(*fi).C() = vcg::Color4b::LightGreen;
		}
	}

	for (MyMeshOcf::FaceIterator fi = m2.face.begin(); fi != m2.face.end(); ++fi) if (!(*fi).IsD()) {
		if ((*fi).V(0)->IsS() || (*fi).V(1)->IsS() || (*fi).V(2)->IsS()) {
			(*fi).C() = vcg::Color4b::LightRed;
		}
	}

	// start with identity as initial transformation
	// in practice you might want to use some kind of prediction here
	Matrix R = Matrix::eye(3);
	Matrix t(3, 1);

	// run point-to-plane ICP (-1 = no outlier threshold)
	cout << endl << "Running ICP (point-to-plane, no outliers)" << endl;

	float t1 = clock();
	IcpPointToPlane icp(M, num, dim);
	double residual = icp.fit(T, num, R, t, -1);
	cout << "ICP process : " << (clock() - t1) / CLOCKS_PER_SEC << endl;

	// results
	cout << endl << "Transformation results:" << endl;
	cout << "R:" << endl << R << endl << endl;
	cout << "t:" << endl << t << endl << endl;
	cout << "Residual:" << residual << endl;

	Matrix Rt = Matrix::inv(R);
	vcg::Matrix44f transM;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			transM[i][j] = R.val[i][j];
		}
		transM[i][3] = t.val[0][i];
	}
	transM[3][0] = 0;
	transM[3][1] = 0;
	transM[3][2] = 0;
	transM[3][3] = 1.0;

	vcg::tri::io::ExporterSTL<MyMeshOcf>::Save(m1, "../result_stl/translated2.stl", true, vcg::tri::io::Mask::IOM_FACECOLOR, 0, 1);
	vcg::tri::UpdatePosition<MyMeshOcf>::Matrix(m2, transM);
	vcg::tri::io::ExporterSTL<MyMeshOcf>::Save(m2, "../result_stl/translated.stl", true, vcg::tri::io::Mask::IOM_FACECOLOR, 0, 1);

	visualizeResult(m2, A, VISUALIZE_ERROR_METRIC);

	// free memory
	free(M);
	free(T);

	// success
	return 0;
}

/* Colorize the distance between the STL and DICOM data.
 * Export colorized mesh for mesh m (STL data).
 * pointset: Dicom Model point set
 * pointToPoint: if true, use point-to-point distance (STL face's centroid to nearest point of DICOM) to show the error
 *               o.w use point-to-plane (distance between the distance of STL face to the nearest point of DICOM to show the error) */
void visualizeResult(MyMeshOcf& m, kdtree::KDTreeArray pointset, bool pointToPoint) {
	// Construct kd tree with dicom points.
	kdtree::KDTree* dicomKDTree = new kdtree::KDTree(pointset);
	std::vector<float>         centroid(3);
	kdtree::KDTreeResultVector neighbor;



	if (!pointToPoint) {
		vcg::tri::UpdateNormal<MyMeshOcf>::PerFaceNormalized(m);
	}


	for (MyMeshOcf::FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD()) {
		centroid[0] = ((*fi).V(0)->P()[0] + (*fi).V(1)->P()[0] + (*fi).V(2)->P()[0]) / 3;
		centroid[0] = ((*fi).V(0)->P()[1] + (*fi).V(1)->P()[1] + (*fi).V(2)->P()[1]) / 3;
		centroid[0] = ((*fi).V(0)->P()[2] + (*fi).V(1)->P()[2] + (*fi).V(2)->P()[2]) / 3;

		// Find the nearest point of the face's centroid.
		dicomKDTree->n_nearest(centroid, 1, neighbor);
		double dicomPx = dicomKDTree->the_data[neighbor[0].idx][0];
		double dicomPy = dicomKDTree->the_data[neighbor[0].idx][1];
		double dicomPz = dicomKDTree->the_data[neighbor[0].idx][2];


		// Distance between stl model's face to the nearest dicom model's point.
		double dist;
		if (pointToPoint) {
			dist = sqrt(pow(dicomPx - centroid[0], 2) + pow(dicomPy - centroid[1], 2) + pow(dicomPx - centroid[2], 2));

		}
		else {
			dist = (*fi).N()[0] * (dicomPx - centroid[0]) +
				(*fi).N()[1] * (dicomPx - centroid[1]) +
				(*fi).N()[2] * (dicomPx - centroid[2]);
		}

		if (dist < -10 * ERROR_SCALE) {
			(*fi).C() = vcg::Color4b::Red;
		}
		else if (dist < -5 * ERROR_SCALE) {
			(*fi).C() = vcg::Color4b::Blue;
		}
		else if (dist < -1 * ERROR_SCALE) {
			(*fi).C() = vcg::Color4b::Green;
		}
		else if (dist < 1 * ERROR_SCALE) {
			(*fi).C() = vcg::Color4b::LightGreen;
		}
		else if (dist < 5 * ERROR_SCALE) {
			(*fi).C() = vcg::Color4b::LightBlue;
		}
		else if (dist < 10 * ERROR_SCALE) {
			(*fi).C() = vcg::Color4b::LightRed;
		}
		else {
			(*fi).C() = vcg::Color4b::LightGray;
		}
	}

	vcg::tri::io::ExporterSTL<MyMeshOcf>::Save(m, "../result_stl/colormap.stl", true, vcg::tri::io::Mask::IOM_FACECOLOR, 0, 1);
}