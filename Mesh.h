#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

using Eigen::MatrixXd;
using namespace std;

struct Vertex
{
	float* coords, * normals; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vertices;
	vector< int > triList; 
	vector< int > edgeList; 
	
	Vertex(int i, float* c) : idx(i), coords(c) {};
	float getCoord(int i) const {return coords[i];} ;
	float* getCoords() const {return coords;} ;
};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	vector< int > triList; //adj triangles (at max 2)
	float length;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) {};

	void computeLength(Vertex* v1, Vertex* v2)
	{
		float* coords1 = v1->coords;
		float* coords2 = v2->coords;
		length = sqrt(pow(coords1[0]-coords2[0], 2) + pow(coords1[1]-coords2[1], 2) + pow(coords1[2]-coords2[2], 2));
	}
};

struct Triangle
{
	int idx; //tris[idx]
	int nonSeedIdx;
	int v1i, v2i, v3i;
	float* center, * normal, * color;

	vector< int > edgeList;
	vector< int > triList; //adj triangles (at max 3)
	vector< float > dihAngList; // dihedral angle measure for each adj triangle resp.	
	float probability;	// for rw segmentation
	int seedRegion;		// for rw segmentation
	bool isItSeed;
	bool isNeighborToSeed;
	
	float sdfValue;		// for sdf segmentation

	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {isItSeed=false; isNeighborToSeed=false;} ;
	int getVi(int i) const {if (i==0) return v1i; else if (i==1) return v2i; else return v3i;} ;
	float getSDF() const {return sdfValue;} ;
	void setSDF(float val) {sdfValue = val;} ;
	void computeCentAndNorm(Vertex* v1, Vertex* v2, Vertex* v3);
	int findCommonEdge(Triangle* neighbor);
};

class Mesh
{
private:
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;
	float** distMatrix;	// holds the geodesic distances
	int** pathMatrix;	// holds the geodesic paths
	vector< float > sdfValues; // holds the sdf value for each triangle

	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2, int t_id);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i, int t_id);

public:
	int groupCount;		// holds the number of Segments
	int groupSize;		// holds the sdf interval size
	float* groupUpperBounds; // holds the sdf upper bound for each segment
	vector< float* > rayStart;
	vector< float* > rayEnd;

	Mesh() {} ;
	~Mesh();
	const Vertex* getVertex(int id) {return verts[id];} ;
	const Triangle* getTriangle(int id) {return tris[id];} ;
	const Edge* getEdge(int id) {return edges[id];} ;
	int getNumOfVerts() {return verts.size();} ;
	int getNumOfEdges() {return edges.size();} ;
	int getNumOfTris() {return tris.size();} ;
	int getPath(int from, int to) {return pathMatrix[from][to];} ;
	float getDistance(int from, int to) {return distMatrix[from][to];} ;
	void loadOff(const char* name);
	void calculateGeodesicsByFW(const char* name);
	void calculateGeodesicsByArray(const char* name);
	void cleanGeodesicData();
	void segmentMeshBySDF(int segmentCount);
	void cleanSDFData();
	void detectSeeds(int segmentCount, Triangle* seedTris[], int selectionType);
	void segmentMeshByRW(int segmentCount, int selectionType);
	float* intersectRayTriangle(float* start, float* dir);
	int getNumOfRays() {return rayStart.size();} ;
	float* getRayStart(unsigned int r) {return rayStart[r];} ;
	float* getRayEnd(unsigned int r) {return rayEnd[r];} ;
};

