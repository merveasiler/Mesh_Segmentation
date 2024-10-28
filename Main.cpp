#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoCube.h>

#include "Mesh.h"
#include "Heap.h"
#include "Painter.h"

int main(int argc, char ** argv)
{
	QWidget * mainwin = SoQt::init(argc, argv, argv[0]);
	SoQtExaminerViewer * viewer = new SoQtExaminerViewer(mainwin);
	SoSeparator * root = new SoSeparator;
	root->ref();

	// get file name
	cout << "Enter the input file name: ";
	string inputFile;
	getline(cin, inputFile);

	// process file
	Painter* painter = new Painter();
	SoSeparator * res;
	Mesh* mesh = new Mesh();
	mesh->loadOff(inputFile.c_str());

	// get the operation type
	cout << "Enter operation type <geodesic: 0, sdf_seg: 1, rw_seg: 2>: ";
	int opType;
	cin >> opType;

	if (opType == 0) {		// DRAW GEODESICS
		// using heap
		struct Graph* graph = createGraph( mesh->getNumOfVerts() );
		for (unsigned int i=0; i < mesh->getNumOfEdges(); i++) {
			const Edge* edge = mesh->getEdge(i);
			addEdge(graph, edge->v1i, edge->v2i, edge->length);
			addEdge(graph, edge->v2i, edge->v1i, edge->length);
		}
		calculateGeodesicByHeap(graph, ("geodesicLengthsByHeap"+inputFile+".txt").c_str());

		// using array
		mesh->calculateGeodesicsByArray(("geodesicLengthsByArray"+inputFile+".txt").c_str());

		// using dynamic programming (Floyd Warshall)
		mesh->calculateGeodesicsByFW(("geodesicLengthsByFW"+inputFile+".txt").c_str());

		res = painter->getShapeSep(mesh, false, false);
		painter->drawTriangulation(res, mesh);
		cout << "From <vertex id> to <vertex id>: ";
		int from, to;
		cin >> from >> to;
		cout << "The distance from the vertex" << from << " to the vertex " << to << " is: " << mesh->getDistance(from, to) << endl;
		painter->drawGeodesic(res, mesh, from, to);

	}
	else {
		cout << "Enter the segment count: ";
		int segmentCount;
		cin >> segmentCount;

		if (opType == 1) {	// SEGMENT BY SHAPE DIAMETER FUNCTION
			mesh->segmentMeshBySDF(segmentCount);
			res = painter->getShapeSep(mesh, true, false);
			cout << "Do you want to see the the inverse of the rays? <yes: y, no: n>: ";
			char raysOnOff;
			cin >> raysOnOff;
			if (raysOnOff == 'y') {
				painter->drawTriangulation(res, mesh);
				for (unsigned r=0; r<mesh->getNumOfRays(); r++)
					painter->drawRays(res, mesh->getRayStart(r), mesh->getRayEnd(r));
			}
		}
		else {			// SEGMENT BY RANDOM WALK
			cout << "What kind of seed selection do you want? <random: 1, fps: 2>: ";
			int selectionType;
			cin >> selectionType;
			mesh->segmentMeshByRW(segmentCount, selectionType);
			res = painter->getShapeSep(mesh, false, true);
		}

	}

	// display on the screen
	root->addChild(res);
	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);	
	viewer->show();
	SoQt::show(mainwin);
	SoQt::mainLoop();
	delete viewer;
	root->unref();

	//mesh->cleanGeodesicData();
	//mesh->cleanSDFData();
	delete painter;
	delete mesh;
	return 0;
} 
