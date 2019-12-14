#pragma once
#include <vcg/complex/complex.h>

/* Opitonal Component Fast of vcg */
class MyVertexOcf;
class MyFaceOcf;
class MyEdge;

struct MyUsedTypesOcf : public vcg::UsedTypes<
	vcg::Use<MyVertexOcf>   ::AsVertexType,
	vcg::Use<MyEdge>     ::AsEdgeType,
	vcg::Use<MyFaceOcf>     ::AsFaceType> {};

class MyVertexOcf : public vcg::Vertex< MyUsedTypesOcf, vcg::vertex::InfoOcf, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Curvatured, vcg::vertex::Color4b, vcg::vertex::BitFlags, vcg::vertex::MarkOcf > {};
class MyFaceOcf : public vcg::Face< MyUsedTypesOcf, vcg::face::InfoOcf, vcg::face::FFAdjOcf, vcg::face::VertexRef, vcg::face::Color4b, vcg::face::BitFlags, vcg::face::Normal3f, vcg::face::MarkOcf > {};
class MyEdge : public vcg::Edge< MyUsedTypesOcf> {};
class MyMeshOcf : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<MyVertexOcf>, vcg::face::vector_ocf<MyFaceOcf> > {};
