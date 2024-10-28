#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoIndexedLineSet.h>

#include "Mesh.h"


class Painter
{
public:
	~Painter();
	SoSeparator* getShapeSep(Mesh* mesh, bool sdfSegmentPainting, bool rwSegmentPainting);
	void drawTriangulation(SoSeparator* res, Mesh* mesh);
	void drawGeodesic(SoSeparator* res, Mesh* mesh, int from, int to);
	void drawNormals(SoSeparator* res, Mesh* mesh);
	void drawRays(SoSeparator* res, float src[3], float dest[3]);
};
