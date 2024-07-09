
#include "alg.hpp"
#include "objutils.hpp"
#include <algorithm>

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

vector<vec2> advancingFrontFromObj(Mesh* mesh) {
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

vector<vector<vec2>> advancingFrontObjFirstStep(Mesh* mesh) {
	vector<vector<vec2>> resultingTriangles{};

	for (MeshObject* meshObject : mesh->objects) {
		//deconstruct each "object", just get all the points
		vector<vec2>* points = new vector<vec2>();
		for (unsigned int index : meshObject->vertexIndices) {
			points->push_back(mesh->vertices.at(index));
		}

		//do the jarvis method for hull
		vector<vec2> convexHull = jarvis(points);

		vector<vec2> triangulatedMeshObject = adf2(*points, convexHull);

		//resultingTriangles.push_back(vector<vec2>());
		//for (auto p : triangulatedMeshObject) {
			//resultingTriangles.push_back(triangulatedMeshObject);
		//}
		resultingTriangles.push_back(triangulatedMeshObject);
	}

	return resultingTriangles;
}

Mesh* advancingFrontObjSecondStep(vector<vector<vec2>> triangles, string meshName) {
	Mesh* finalMesh = new Mesh(meshName);

	//reconstruct the object
	for (vector<vec2> objectPointsVector : triangles) {
		cout << "new object" << endl;
		//cout << "a" << endl;
		MeshObject* currentObject = new MeshObject();
		MeshFace* currentFace = new MeshFace();
		int pointCount = 0; // when this hits multiples of 3, it's a face.
		for (vec2 point : objectPointsVector) {
			
			// vertices can be repeated in the mesh if they are from different objects
			int foundIndex = -1;
			for (int testIndex : currentObject->vertexIndices) {
				if (point == finalMesh->vertices[testIndex]) {
					foundIndex = testIndex;
					break;
				}
			}
			// if point doesn't already exist, add it to the end of mesh and object
			if (foundIndex == -1) {
				finalMesh->vertices.push_back(point);
				foundIndex = finalMesh->vertices.size() - 1;
				currentObject->vertexIndices.push_back(foundIndex);
				//cout << "added" << endl;
			}

			currentFace->vertexIndices.push_back(foundIndex);

			// create another face if we processed 3 points already
			pointCount += 1;
			if (pointCount % 3 == 0) {
				currentObject->faces.push_back(currentFace);
				cout << currentFace->vertexIndices[0]+1 << " " << currentFace->vertexIndices[1]+1 << " " << currentFace->vertexIndices[2]+1 << endl;
				currentFace = new MeshFace();
			}
		}

		finalMesh->objects.push_back(currentObject);
	}

	return finalMesh;
}

//Mesh* triangulatedToMesh(vector<vec2> points) {
//	Mesh* finalMesh = new Mesh();
//
//	for ()
//}