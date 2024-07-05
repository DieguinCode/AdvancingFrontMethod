
#include "alg.hpp"
#include "objutils.hpp"


using namespace std;

Mesh* jarvisObj(Mesh* mesh) {
	Mesh* convexedMesh = new Mesh();

	unsigned int vertexCounter = 0; //global vertex counter
	for (MeshObject* meshObject : mesh->objects) {

		//deconstruct each "object", just get all the points
		vector<vec2>* points = new vector<vec2>();
		for (unsigned int index : meshObject->vertexIndices) {
			points->push_back(mesh->vertices.at(index));
		}

		//do the jarvis method
		vector<vec2> convexHull = jarvis(points);

		//reconstruct a single face with all the points
		MeshObject* newObject = new MeshObject();
		newObject->name = meshObject->name;
		MeshFace* newFace = new MeshFace();
		newObject->faces.push_back(newFace);

		for (vec2 point : convexHull) {
			convexedMesh->vertices.push_back(point);
			newObject->vertexIndices.push_back(vertexCounter);
			newFace->vertexIndices.push_back(vertexCounter);
			vertexCounter += 1; //global vertex counter
		}

		//save each individual object
		convexedMesh->objects.push_back(newObject);
	}

	return convexedMesh;
}

vector<vec2> advancingFrontObj(Mesh* mesh) {
	vector<vec2> resultingTriangles{};

	for (MeshObject* meshObject : mesh->objects) {
		//deconstruct each "object", just get all the points
		vector<vec2>* points = new vector<vec2>();
		for (unsigned int index : meshObject->vertexIndices) {
			points->push_back(mesh->vertices.at(index));
		}

		//do the jarvis method for hull
		vector<vec2> convexHull = jarvis(points);

		vector<vec2> triangulatedMeshObject = adf2(*points, convexHull);

		for (auto p : triangulatedMeshObject) {
			resultingTriangles.push_back(p);
		}
	}

	return resultingTriangles;
}