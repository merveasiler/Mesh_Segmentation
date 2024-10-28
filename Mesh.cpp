#include "Mesh.h"

void Mesh::detectSeeds(int segmentCount, Triangle* seedTris[], int selectionType) {

	if (selectionType == 1) { // random seed selection
		for (int num=0; num<segmentCount; num++) {
			int id = rand() % tris.size();
			bool valid = true;
			for (int s=0; s<num; s++)
				if (seedTris[s]->idx == id) {
					num--;
					valid = false;
					break;
				}
			if (valid)
				seedTris[num] = tris[id];
		}
		return;
	}

	// find the centroid of the mesh
	float centroid[3] = {0,0,0};
	for (unsigned int i=0; i < tris.size(); i++) {
		centroid[0] += tris[i]->center[0];
		centroid[1] += tris[i]->center[1];
		centroid[2] += tris[i]->center[2];
	}
	
	centroid[0] /= tris.size();
	centroid[1] /= tris.size();
	centroid[2] /= tris.size(); 

	// find the triangle closest to centroid
	Triangle* central_t;
	float shortestDist = -1, dist;
	for (unsigned int i=0; i < tris.size(); i++) {
		dist = sqrt( pow((tris[i]->center[0] - centroid[0]), 2) + 
			     pow((tris[i]->center[1] - centroid[1]), 2) + 
			     pow((tris[i]->center[2] - centroid[2]), 2) );
		if (dist < shortestDist || shortestDist < 0) {
			shortestDist = dist;
			central_t = tris[i];
		}
	}

	// find the first seed; the farthest to central_t
	calculateGeodesicsByFW("temp.txt");
	float longestDist = 0;
	for (unsigned int i=0; i < verts.size(); i++) {
		dist = distMatrix[central_t->v1i][i];
		if (dist > longestDist) {
			longestDist = dist;
			seedTris[0] = tris[ verts[i]->triList[0] ];
		}
	}

	// find the other seeds
	for (int num=1; num<segmentCount; num++) {	// for each segment
		
		longestDist = 0;
		for (unsigned int i=0; i < verts.size(); i++) {	// check each vertex
			bool skip = false;
			for (int s=0; s<num; s++) {	// dis we choose vertex of some seed?
				if (seedTris[s]->v1i == i || seedTris[s]->v2i == i || seedTris[s]->v3i == i)
					{ skip = true; break; }
			}
			if (skip) continue;

			shortestDist = -1;
			for (int s=0; s<num; s++) {	// previous seeds
				dist = distMatrix[i][seedTris[s]->v1i] ;
				if (dist < shortestDist || shortestDist < 0)
					shortestDist = dist;
			}
			
			if (shortestDist > longestDist) {
				longestDist = shortestDist;
				seedTris[num] = tris[ verts[i]->triList[0] ];
			}
		}
	}

}

void Mesh::segmentMeshByRW(int segmentCount, int selectionType) {
	/************************ initializations ***********************/
	groupCount = segmentCount;
	Triangle* seedTris[groupCount];
	detectSeeds(groupCount, seedTris, selectionType);
	// kuyruk, sag arka ayak, sol arka ayak, sag on ayak, sol on ayak, boyun, govde
	//Triangle* seedTris[groupCount] = {tris[999], tris[300], tris[400], tris[900], tris[990], tris[700], tris[175] };
	for (unsigned int i=0; i<groupCount; i++) {
		Triangle* seed = seedTris[i];
		seed->probability = 1;
		seed->seedRegion= i;
		seed->isItSeed = true;
		for (unsigned int j=0; j < seed->triList.size(); j++) {
			Triangle* seedNeighbor = tris[seed->triList[j]];
			seedNeighbor->isNeighborToSeed = true;
			seedNeighbor->seedRegion = i;
			seedNeighbor->probability = 1;
		}
	}

	// extract non-seed triangles
	vector< Triangle* > nonSeedTris;
	for (unsigned int idx=0, index=0; idx<tris.size(); idx++) {
		if (tris[idx]->isItSeed)
			continue;
		nonSeedTris.push_back(tris[idx]);
		tris[idx]->nonSeedIdx = index;
		index++;
	}
	/**************************************************************/

	// calculate dihedral angle measures
	for (unsigned int i=0; i<tris.size(); i++) {
		Triangle* t = tris[i];
		float totalDihAng = 0;

		for (unsigned int j=0; j < t->triList.size(); j++) {
			Triangle* adj = tris[ t->triList[j] ];
			float cosDihedral = -1*(adj->normal[0]*t->normal[0] + adj->normal[1]*t->normal[1] + adj->normal[2]*t->normal[2]);

			float nuCoeff = 1; // concave
			float centralAxis[3] = {adj->center[0]-t->center[0], adj->center[1]-t->center[1], adj->center[2]-t->center[2]};
			float length = sqrt(centralAxis[0]*centralAxis[0] + centralAxis[1]*centralAxis[1] + centralAxis[2]*centralAxis[2]);
			centralAxis[0] = centralAxis[0]/length; centralAxis[1] = centralAxis[1]/length; centralAxis[2] = centralAxis[2]/length;
			if (centralAxis[0]*adj->normal[0] + centralAxis[1]*adj->normal[1] + centralAxis[2]*adj->normal[2] > 0) // convex
				nuCoeff = 0.1;

			float d_a = nuCoeff * (1-cosDihedral);
			t->dihAngList.push_back( d_a );
			totalDihAng += d_a;
		}

		// normalize
		for (unsigned int j=0; j < t->triList.size(); j++)
			t->dihAngList[j] = t->dihAngList[j] / totalDihAng;
	}

	// solve the linear system for each seed
	for (unsigned int i=0; i<groupCount; i++) { // walk through each seed
		Triangle* seed = seedTris[i];
		Eigen::SparseMatrix<double> A(nonSeedTris.size(), nonSeedTris.size());
		Eigen::VectorXd B(nonSeedTris.size());

		// prepare A and B matrices
		for (unsigned int l=0; l<nonSeedTris.size(); l++) { // walk through non-seed
			Triangle* t = nonSeedTris[l];
			A.coeffRef(l,l) = 1;
			B(l) = 0;
			vector<float> probList;
			float totalProb = 0;

			for (unsigned int n=0; n < t->triList.size(); n++) { // walk through neighbors-1
				Triangle* neighbor = tris[t->triList[n]];
				Edge* commonEdge = edges[t->findCommonEdge(neighbor)];
				float probability = (commonEdge->length)*exp( -1*(t->dihAngList[n]) );
				probList.push_back(probability);
				totalProb += probability;
			}
			for (unsigned int n=0; n < t->triList.size(); n++) { // walk through neighbors-2
				Triangle* neighbor = tris[t->triList[n]];
				if (neighbor->isItSeed) {
					if (neighbor->idx == seed->idx)
						B(l) += (probList[n]/totalProb);
				}
				else
					A.coeffRef(l, neighbor->nonSeedIdx) = -(probList[n]/totalProb);
			}
			probList.clear();
		}

		// solve the system
		A.makeCompressed();
		Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > qr(A);
		Eigen::VectorXd P = qr.solve(B);

		// save the value
		for (unsigned int l=0; l<nonSeedTris.size(); l++) {
			Triangle* t = nonSeedTris[l];
			if ( P(l) > t->probability ) {
				t->probability = P(l);
				t->seedRegion = i;
			}
		}
	}

	nonSeedTris.clear();
}

void Mesh::segmentMeshBySDF(int segmentCount) {

	/****************** SDF Method *****************/

	float pi = 3.14159265;
	int angleCount = 5;
	int maxAngle = 45;
	float angleInterval = maxAngle / angleCount;
	float tangents[angleCount];
	float weights[angleCount];
	for (int i=1; i<=angleCount; i++) {
		tangents[i] = tan(angleInterval*i*pi/180);
		weights[i] = cos(angleInterval*i*pi/180);
	}
	float points[6][3];
	float target[3]; // indicates target point giving the direction of the ray starting at center
	//float dir[28][3]; // each row indicates the direction of the ray
	float dir[3];

	for (unsigned int i=0; i<tris.size(); i++) {

		float _SDF = 0;
		int rayCount = 0;

		Triangle* t = tris[i];
		float* center = t->center;
		float innerNormal[3] = {-(t->normal[0]), -(t->normal[1]), -(t->normal[2])};
		float p[3] = {center[0]+innerNormal[0], center[1]+innerNormal[1], center[2]+innerNormal[2]}; 

		// define the points that the ray in the plane will go to from center
		for (unsigned int j=0; j<3; j++) {
			points[0][j] = verts[t->getVi(0)]->coords[j];
			points[1][j] = verts[t->getVi(1)]->coords[j];
			points[2][j] = verts[t->getVi(2)]->coords[j];
			points[3][j] = (points[0][j] + points[1][j]) / 2;
			points[4][j] = (points[0][j] + points[2][j]) / 2;
			points[5][j] = (points[1][j] + points[2][j]) / 2;
		}

		// define the ray directions
		for (unsigned int j=0; j<6; j++) {
			float vec[3] = {points[j][0] - center[0], points[j][1] - center[1], points[j][2] - center[2]};
			float height = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
			vec[0] = vec[0] / height, vec[1] = vec[1] / height, vec[2] = vec[2] / height;
			for (unsigned int k=0; k<angleCount; k++) {
				// compute the ray direction
				target[0] = p[0] + tangents[k]*vec[0];
				target[1] = p[1] + tangents[k]*vec[1];
				target[2] = p[2] + tangents[k]*vec[2];
				dir[0] = target[0]-center[0], dir[1] = target[1]-center[1], dir[2] = target[2]-center[2];
				float length = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
				dir[0] = dir[0] / length, dir[1] = dir[1] / length, dir[2] = dir[2] / length;
			
				// save rays to draw later
				float* outTarget = new float[3];
				outTarget[0] = center[0]-5*dir[0]; outTarget[1] = center[1]-5*dir[1], outTarget[2] = center[2]-5*dir[2];
				rayStart.push_back(center);
				rayEnd.push_back(outTarget);

				// compute the ray triangle intersection
				float* intersection = intersectRayTriangle(center, dir);
				if (intersection[3] >= 0) {	// an intersection exists
					float* n = tris[intersection[3]]->normal;	// compute SDF value
					if ( acos(n[0]*dir[0] + n[1]*dir[1] + n[2]*dir[2]) <= pi/2 ) { // else: ignore
						_SDF += intersection[4]*weights[k];	
						rayCount++;
					}
				}
				delete[] intersection;
			}
		}

		// compute the triangle intersection for the ray in direction of innerNormal
		dir[0] = innerNormal[0]; dir[1] = innerNormal[1]; dir[2] = innerNormal[2];
		float* intersection = intersectRayTriangle(center, dir);
		if (intersection[3] >= 0) {	// an intersection exists
			float* n = tris[intersection[3]]->normal;	// compute SDF value
			if ( acos(n[0]*dir[0] + n[1]*dir[1] + n[2]*dir[2]) <= pi/2 ) { // else: ignore
				_SDF += intersection[4];	
				rayCount++;
			}
		}
		delete[] intersection;

		if (rayCount > 0)
			sdfValues.push_back(_SDF/rayCount);
		else	sdfValues.push_back(0);
		t->sdfValue = sdfValues[i];
	}

	// construct histogram
	groupCount = segmentCount;
	std::sort(sdfValues.begin(), sdfValues.begin()+sdfValues.size());
	float diff = sdfValues[sdfValues.size()-1] - sdfValues[0];
	groupSize = (diff / groupCount) + 1;
	groupUpperBounds = new float[groupCount];
	for (unsigned int i=0; i<groupCount; i++)
		groupUpperBounds[i] = sdfValues[0]+(i+1)*groupSize;
}

void Mesh::cleanSDFData() {
	delete[] groupUpperBounds;
	sdfValues.clear();
}

void Mesh::calculateGeodesicsByArray(const char* outputFile) {

	int i, j, k;
	int nVerts = verts.size(), nEdges = edges.size();
	ofstream myFile (outputFile);
	clock_t timeMeasure;
	timeMeasure = clock();

	for (int src=0; src<nVerts; src++) {
		float distMatrix[nVerts];
		float pathMatrix[nVerts];

		for (i=0; i<nVerts; i++)
			if (i == src) {
				distMatrix[i] = 0;
				pathMatrix[i] = i;
			}
			else {
				distMatrix[i] = -1;
				pathMatrix[i] = -1;
			}
	}

	myFile.close();

}

void Mesh::calculateGeodesicsByFW(const char* outputFile) {

	int i, j, k;
	int nVerts = verts.size(), nEdges = edges.size();
	ofstream myFile (outputFile);
	clock_t timeMeasure;
	timeMeasure = clock();

	distMatrix = new float*[nVerts];
	pathMatrix = new int*[nVerts];
	for (i=0; i<nVerts; i++) {
		distMatrix[i] = new float[nVerts];
		pathMatrix[i] = new int[nVerts];
	}

	for (i=0; i<nVerts; i++)
		for (j=0; j<nVerts; j++)
			if (i == j) {
				distMatrix[i][j] = 0;
				pathMatrix[i][j] = i;
			}
			else {
				distMatrix[i][j] = -1;
				pathMatrix[i][j] = -1;
			}

	for (i=0; i<nEdges; i++) {
		Edge* edge = edges[i];
		distMatrix[edge->v1i][edge->v2i] = edge->length;
		distMatrix[edge->v2i][edge->v1i] = edge->length;
		pathMatrix[edge->v1i][edge->v2i] = edge->v1i;
		pathMatrix[edge->v2i][edge->v1i] = edge->v2i;
	}

	for (k=0; k<nVerts; k++)
		for (i=0; i<nVerts; i++)
			for (j=0; j<nVerts; j++) {
				if (distMatrix[i][k] < 0 || distMatrix[k][j] < 0)
					continue;
				if (distMatrix[i][j] < 0 || distMatrix[i][j] > distMatrix[i][k] + distMatrix[k][j]) {
					distMatrix[i][j] = distMatrix[i][k] + distMatrix[k][j];
					pathMatrix[i][j] = pathMatrix[k][j];
				}
			}

	for (i=0; i<nVerts; i++) {
		for (j=0; j<nVerts; j++)
			myFile << distMatrix[i][j] << " ";
		myFile << endl;
	}
	myFile.close();
	timeMeasure = clock() - timeMeasure;
	cout << "It took " << ((float)timeMeasure)/CLOCKS_PER_SEC << "seconds to calculate all geodesics by Floyd Warshall.\n";
}

void Mesh::cleanGeodesicData() {
	for (unsigned int i=0; i<verts.size(); i++) {
		delete[] distMatrix[i];
		delete[] pathMatrix[i];
	}
	
	delete[] distMatrix;
	delete[] pathMatrix;
}

void Mesh::loadOff(const char* name)
{
	int nVerts, nTris, n, id, i = 0;
	float x, y, z;
	string line;

	ifstream myFile (name);
	if (myFile.is_open())
	{
		getline (myFile,line);	// read "OFF"
		myFile >> nVerts >> nTris >> n;

		while (i++ < nVerts)
		{
			myFile >> x >> y >> z;
			addVertex(x, y, z);
		}

		i=0;
		while (i++ <nTris)
		{
			myFile >> id >> x >> y >> z;
			addTriangle((int) x, (int) y, (int) z);
		}
		myFile.close();
	}
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back( new Vertex(idx, c) );
}


void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	Triangle* t = new Triangle(idx, v1, v2, v3);
	t->computeCentAndNorm(verts[v1], verts[v2], verts[v3]);
	tris.push_back(t);
	t->probability = -1;	// default value for rw segmentation
	t->seedRegion = 0;	// default value for rw segmentation
	t->color = new float[3];

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	if (! makeVertsNeighbor(v1, v2, idx) )
		addEdge(v1, v2, idx);

	if (! makeVertsNeighbor(v1, v3, idx) )
		addEdge(v1, v3, idx);

	if (! makeVertsNeighbor(v2, v3, idx) )
		addEdge(v2, v3, idx);

}

bool Mesh::makeVertsNeighbor(int v1i, int v2i, int t_id)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->edgeList.size(); i++) {
		Edge* edge = edges[verts[v1i]->edgeList[i]];
		if (edge->v1i == v2i || edge->v2i == v2i) {

			Triangle* tNew = tris[t_id];
			tNew->edgeList.push_back(edge->idx);

			for (int j=0; j < edge->triList.size(); j++) {
				Triangle* tOld = tris[edge->triList[j]];
				tOld->triList.push_back(tNew->idx);
				tNew->triList.push_back(tOld->idx);
			}

			edge->triList.push_back(t_id);
			return true;
/*
			Triangle* t1 = tris[edge->triList[0]];
			Triangle* t2 = tris[t_id];

			edge->triList.push_back(t_id);
			t2->edgeList.push_back(edge->idx);

			t1->triList.push_back(t2->idx);
			t2->triList.push_back(t1->idx);
			return true;
*/
		}
	}


	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addEdge(int v1, int v2, int t_id)
{
	int idx = edges.size();
	edges.push_back( new Edge(idx, v1, v2) );
	tris[t_id]->edgeList.push_back(idx);
	edges[idx]->triList.push_back(t_id);
	edges[idx]->computeLength(verts[v1], verts[v2]);

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

Mesh::~Mesh() {
	int i;
	for (i=0; i<tris.size(); i++) {
		delete[] tris[i]->center;
		delete[] tris[i]->normal;
		delete[] tris[i]->color;
		(tris[i]->triList).clear();
		(tris[i]->edgeList).clear();
		(tris[i]->dihAngList).clear();
		delete tris[i];
	}

	for (i=0; i<edges.size(); i++) {
		(edges[i]->triList).clear();
		delete edges[i];
	}

	for (i=0; i<verts.size(); i++) {
		(verts[i]->vertList).clear();
		(verts[i]->triList).clear();
		(verts[i]->edgeList).clear();
		delete[] verts[i]->coords;
		delete verts[i];
	}

	for (i=0; i<rayEnd.size(); i++)
		delete[] rayEnd[i];
	rayEnd.clear();
}

float computeDet(float* column1, float* column2, float* column3) {
	float part1 = column1[0] * (column2[1]*column3[2] - column2[2]*column3[1]);
	float part2 = -column2[0] * (column1[1]*column3[2] - column1[2]*column3[1]);
	float part3 = column3[0] * (column1[1]*column2[2] - column1[2]*column2[1]);
	return part1+part2+part3;
}

float* Mesh::intersectRayTriangle(float* start, float* dir) {
	
	float* intersection = new float[5]; // the intersection point [0,1,2], triangle id [3] and tVal [4], resp.
	intersection[0] = 0, intersection[1] = 0, intersection[2] = 0, intersection[3] = -1, intersection[4] = -1;

	for (unsigned int i=0; i<tris.size(); i++) {

		Triangle* t = tris[i];
		float* a = verts[t->getVi(0)]->coords;
		float* b = verts[t->getVi(1)]->coords;
		float* c = verts[t->getVi(2)]->coords;

		float resultColumn[3] = {a[0]-start[0], a[1]-start[1], a[2]-start[2]};
		float column1[3] = {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
		float column2[3] = {a[0]-c[0], a[1]-c[1], a[2]-c[2]};

		float detA = computeDet(column1, column2, dir);
		float detBeta = computeDet(resultColumn, column2, dir);
		float detGama = computeDet(column1, resultColumn, dir);
		float det_t = computeDet(column1, column2, resultColumn);

		float beta = detBeta / detA;
		float gama = detGama / detA;
		float alpha = 1 - (beta+gama);
		float tVal = det_t / detA;
		float epsilon = 0.000001;

		if (alpha >= epsilon && alpha <= 1 && beta >= epsilon && beta <= 1 && gama >= epsilon && gama <= 1 && tVal>0) {
			intersection[0] = start[0] + dir[0]*tVal;
			intersection[1] = start[1] + dir[1]*tVal;
			intersection[2] = start[2] + dir[2]*tVal;
			intersection[3] = i;
			intersection[4] = tVal;
			return intersection;
		}
	}

	return intersection;
}

void Triangle::computeCentAndNorm(Vertex* v1, Vertex* v2, Vertex* v3) {
	// compute center
	float x = (v1->coords[0] + v2->coords[0])/2;
	float y = (v1->coords[1] + v2->coords[1])/2;
	float z = (v1->coords[2] + v2->coords[2])/2;

	center = new float[3];
	center[0] = v3->coords[0] + (x - v3->coords[0])*2/3.0;
	center[1] = v3->coords[1] + (y - v3->coords[1])*2/3.0;
	center[2] = v3->coords[2] + (z - v3->coords[2])*2/3.0;

	// compute normal
	float x1 = v2->coords[0] - v1->coords[0];
	float y1 = v2->coords[1] - v1->coords[1];
	float z1 = v2->coords[2] - v1->coords[2];
	float x2 = v3->coords[0] - v1->coords[0];
	float y2 = v3->coords[1] - v1->coords[1];
	float z2 = v3->coords[2] - v1->coords[2];

	x = y1*z2 - y2*z1;
	y = z1*x2 - z2*x1;
	z = x1*y2 - x2*y1;
	float length = sqrt(x*x + y*y + z*z);

	normal = new float[3];
	normal[0] = x/length;
	normal[1] = y/length;
	normal[2] = z/length;
}

int Triangle::findCommonEdge(Triangle* neighbor) {
	for (int i=0; i < edgeList.size(); i++) {
		int edgeId = edgeList[i];
		for (int j=0; j < neighbor->edgeList.size(); j++)
			if (edgeId == neighbor->edgeList[j])
				return edgeId;
	}
}

