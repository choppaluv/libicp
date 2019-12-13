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
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/math/matrix44.h>
#include "icpPointToPlane.h"
#include "MeshModel.h"

using namespace std;

int main (int argc, char** argv) {

  // define a 3 dim problem with 10000 model points
  // and 10000 template points:
  // This configuration can be changed in config.h file.
  int32_t dim = DIM;
  int32_t num = NUM_SAMPLE_POINTS;

  // allocate model and template memory
  double* M = (double*)calloc(3*num,sizeof(double));
  double* T = (double*)calloc(3*num,sizeof(double));

  MyMesh m1, m2;
  MyMesh::FaceIterator f1;
  const char *filename1 = "../source_stl/LowerJaw.stl";
  const char *filename2 = "../source_stl/DicomModel.stl";

  // import STL file
  int dummy = 0;
  if (vcg::tri::io::ImporterSTL<MyMesh>::Open(m1, filename1, dummy, 0) != 0) {
	  std::printf("Error reading file  %s\n", filename1);
	  exit(0);
  }

  if (vcg::tri::io::ImporterSTL<MyMesh>::Open(m2, filename2, dummy, 0) != 0) {
	  std::printf("Error reading file  %s\n", filename2);
	  exit(0);
  }

  double* A = (double*)calloc(3 * m1.VN(), sizeof(double));
  double* B = (double*)calloc(3 * m2.VN(), sizeof(double));

  int i = 0;
  for (MyMesh::VertexIterator vi = m1.vert.begin();vi != m1.vert.end();++vi) {
	  A[i*3+0] = (*vi).P()[0];
	  A[i*3+1] = (*vi).P()[1];
	  A[i*3+2] = (*vi).P()[2];
	  i++;
  }
  int numVertex1 = i;

  i = 0;
  for (MyMesh::VertexIterator vi = m2.vert.begin();vi != m2.vert.end();++vi) {
	  B[i*3+0] = (*vi).P()[0];
	  B[i*3+1] = (*vi).P()[1];
	  B[i*3+2] = (*vi).P()[2];
	  i++;
  }
  int numVertex2 = i;

  std::random_device randDev;
  std::mt19937 generator(randDev());
  std::uniform_int_distribution<int> distrA(0, numVertex1 - 1);
  std::uniform_int_distribution<int> distrB(0, numVertex2 - 1);
  // int sampleSize = std::min(numVertex1, numVertex2);

  cout << "stl size: " << numVertex1 << endl;
  cout << "dicom size: " << numVertex2 << endl;
  cout << "sample size: " << num << endl;
  
  // set model and template points
  for (int i = 0; i < num; i++) {
	  int idx1 = distrA(generator);
	  int idx2 = distrB(generator);

	  M[i*3+0] = A[idx1*3+0];
	  M[i*3+1] = A[idx1*3+1];
	  M[i*3+2] = A[idx1*3+2];
	  T[i*3+0] = B[idx2*3+0];
	  T[i*3+1] = B[idx2*3+1];
	  T[i*3+2] = B[idx2*3+2];
  }


  // start with identity as initial transformation
  // in practice you might want to use some kind of prediction here
  Matrix R = Matrix::eye(3);
  Matrix t(3,1);

  // run point-to-plane ICP (-1 = no outlier threshold)
  cout << endl << "Running ICP (point-to-plane, no outliers)" << endl;
  IcpPointToPlane icp(M,num,dim);
  double residual = icp.fit(T,num,R,t,-1);

  // results
  cout << endl << "Transformation results:" << endl;
  cout << "R:" << endl << R << endl << endl;
  cout << "t:" << endl << t << endl << endl;
  cout << "Residual:"<<residual;

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


  vcg::tri::UpdatePosition<MyMesh>::Matrix(m2, transM);
  vcg::tri::io::ExporterSTL<MyMesh>::Save(m2, "../result_stl/translated.stl", true, vcg::tri::io::Mask::IOM_FACECOLOR, 0, 1);


  // free memory
  free(M);
  free(T);

  // success
  return 0;
}

