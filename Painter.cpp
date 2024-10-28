#include "Painter.h"

SoSeparator* Painter::getShapeSep(Mesh* mesh, bool sdfSegmentPainting, bool rwSegmentPainting)
{
	SoSeparator* res = new SoSeparator();
	bool segmentPainting = (sdfSegmentPainting || rwSegmentPainting);
	int colors[8][3] = {{1,0,0}, {1,1,0}, {1,0,1}, {0,1,0}, {0,1,1}, {0,0,1}, {1,1,1}, {0,0,0}};

	//color
	SoMaterial* mat = new SoMaterial();

	if (sdfSegmentPainting)	// sdf segmentation colors
		for (unsigned int i=0; i < mesh->getNumOfTris(); i++) {
			float sdfValue = mesh->getTriangle(i)->getSDF();
			// determine color
			int* color;
			for (unsigned j=0; j<mesh->groupCount; j++)
				if (sdfValue < mesh->groupUpperBounds[j]) {
					color = colors[j%8];
					break;
				}
			mat->diffuseColor.set1Value(i, color[0], color[1], color[2]);
		}
	if (rwSegmentPainting)	// rw segmentation colors
		for (unsigned int i=0; i < mesh->getNumOfTris(); i++) {
			int rwValue = mesh->getTriangle(i)->seedRegion;
			int* color = colors[rwValue%8];
			mat->diffuseColor.set1Value(i, color[0], color[1], color[2]);
		}

	res->addChild(mat);

	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints); //Gouraud shading

	if (segmentPainting)
	{
		SoMaterialBinding* materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
		materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
		res->addChild(materialBinding);
	}

	//shape
	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->getNumOfVerts(); c++)
		coords->point.set1Value(c, mesh->getVertex(c)->getCoord(0), mesh->getVertex(c)->getCoord(1), mesh->getVertex(c)->getCoord(2));

	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->getNumOfTris(); c++)
	{
		faceSet->coordIndex.set1Value(c*4, mesh->getTriangle(c)->getVi(0));
		faceSet->coordIndex.set1Value(c*4 + 1, mesh->getTriangle(c)->getVi(1));
		faceSet->coordIndex.set1Value(c*4 + 2, mesh->getTriangle(c)->getVi(2));
		faceSet->coordIndex.set1Value(c*4 + 3, -1);

		if (segmentPainting)
		{
			faceSet->materialIndex.set1Value(0 + 4*c, c);
			faceSet->materialIndex.set1Value(1 + 4*c, c);
			faceSet->materialIndex.set1Value(2 + 4*c, c);
		}
	}

	res->addChild(coords);
	res->addChild(faceSet);

	return res;
}

void Painter::drawTriangulation(SoSeparator* res, Mesh* mesh) {
	// draw triangles
	SoSeparator* thickEdgeSep = new SoSeparator;
	//material
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1,000000;
	thickEdgeSep->addChild(sty);

	//shape
	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;

	//assumes no edge in sedges is removed
	for (unsigned int se = 0; se < mesh->getNumOfEdges(); se++)
	{
		SbVec3f end1 = mesh->getVertex( mesh->getEdge(se)->v1i )->getCoords(),// + SbVec3f(deltaX, 0.0f, 0.0f),
		end2 = mesh->getVertex( mesh->getEdge(se)->v2i )->getCoords();// + SbVec3f(deltaX, 0.0f, 0.0f);
		co->point.set1Value(2*se, end1);
		co->point.set1Value(2*se + 1, end2);
	}

	for (unsigned int ci = 0; ci < mesh->getNumOfEdges(); ci++)
	{
		ils->coordIndex.set1Value(3*ci, 2*ci);
		ils->coordIndex.set1Value(3*ci + 1, 2*ci + 1);
		ils->coordIndex.set1Value(3*ci + 2, -1); //end this edge with -1
	}
	
	thickEdgeSep->addChild(co);
	thickEdgeSep->addChild(ils);
	res->addChild( thickEdgeSep );
	
}

void Painter::drawGeodesic(SoSeparator* res, Mesh* mesh, int from, int to) {

	int k;
	SoSeparator* geodesicSep = new SoSeparator();;
	SoDrawStyle* sty = new SoDrawStyle;
	sty->lineWidth = 1,000000;
	geodesicSep->addChild(sty);

	while(from != to) {
		k = mesh->getPath(from, to);
		const Vertex* v1 = mesh->getVertex(k);
		const Vertex* v2 = mesh->getVertex(to);
		//draw(k,j);
		SoVertexProperty* vprop = new SoVertexProperty();
		vprop->vertex.set1Value(0, v1->coords[0], v1->coords[1], v1->coords[2]);  // Set first vertex to be k
		vprop->vertex.set1Value(1, v2->coords[0], v2->coords[1], v2->coords[2]);  // Set second vertex to be j

		vprop->orderedRGBA = 0xff0000ff;
		SoLineSet* line = new SoLineSet();
		line->vertexProperty = vprop;
		geodesicSep->addChild( line );
		
		to = k;
	}

	res->addChild( geodesicSep );
}

void Painter::drawRays(SoSeparator* res, float src[3], float dest[3]) {

	SoVertexProperty* vprop = new SoVertexProperty();
	vprop->vertex.set1Value(0, src[0], src[1], src[2]);
	vprop->vertex.set1Value(1, dest[0], dest[1], dest[2]);

	vprop->orderedRGBA = 0x00ff0f0f;
	SoLineSet* line = new SoLineSet();
	line->vertexProperty = vprop;
	res->addChild( line );

}

void Painter::drawNormals(SoSeparator* res, Mesh* mesh) {

	for (unsigned int i=0; i<mesh->getNumOfTris(); i++) {
		const Triangle* tris = mesh->getTriangle(i);
		SoVertexProperty* vprop = new SoVertexProperty();
		vprop->vertex.set1Value(0, tris->center[0], tris->center[1], tris->center[2]);
		vprop->vertex.set1Value(1, tris->center[0]+tris->normal[0], 
					tris->center[1]+tris->normal[1], 
					tris->center[2]+tris->normal[2]);

		vprop->orderedRGBA = 0x00ff00ff;
		SoLineSet* line = new SoLineSet();
		line->vertexProperty = vprop;
		res->addChild( line );
	}

}

Painter::~Painter() {
}

