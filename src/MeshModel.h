#pragma once
#include <vcg/complex/complex.h>

class MyVertex;
class MyFace;
class MyEdge;

struct MyUsedTypes : public vcg::UsedTypes<
	vcg::Use<MyVertex>   ::AsVertexType,
	vcg::Use<MyEdge>     ::AsEdgeType,
	vcg::Use<MyFace>     ::AsFaceType> {};

class MyVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Curvatured, vcg::vertex::Color4b, vcg::vertex::BitFlags > {};
class MyFace : public vcg::Face< MyUsedTypes, vcg::face::FFAdj, vcg::face::VertexRef, vcg::face::Color4b, vcg::face::BitFlags, vcg::face::Normal3f > {};
class MyEdge : public vcg::Edge< MyUsedTypes> {};
class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>  > {};