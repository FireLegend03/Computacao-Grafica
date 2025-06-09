#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>  
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

// Default window settings
int width = 800;
int height = 800;

// Default camera settings
float posX = 5.0;
float posY = 5.0;
float posZ = 5.0;

float lookAtX = 0.0;
float lookAtY = 0.0;
float lookAtZ = 0.0;

float upX = 0.0;
float upY = 1.0;
float upZ = 0.0;

float fov = 45.0;
float near = 1.0;
float far = 1000;

float alfa = 0.0f;
float beta = 0.0f;

// Model vars declaration
std::vector<std::string> modelList;
int modelCounter = 0;

// Paths
std::string filePath = "..\\..\\xml\\sistema_solar.xml";
std::string modelsPath = "..\\..\\models\\";

// Configs
bool mousePressed = false;
float movementSpeed = 0.5;

bool flagAxis = true;
bool flagReload = false;
int glMode = GL_LINE;
int glFaces = GL_FRONT;

struct transformation {

	std::string type;

	float angle;

	float x;
	float y;
	float z;

	transformation(std::string type, float angle, float x, float y, float z)
		: type(type), angle(angle), x(x), y(y), z(z) {
	}
};

std::multimap<string, std::vector<transformation>> modelTrans;

void ExplorerModeCam(float radius) {

	posX = radius * cos(beta) * cos(alfa);
	posY = radius * sin(beta);
	posZ = radius * cos(beta) * sin(alfa);
}

void drawAxis(void) {

	glBegin(GL_LINES);

	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(50.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 50.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 50.0f);

	glColor3f(1.0f, 1.0f, 1.0f);

	glEnd();
}

void drawModels(void) {

	int c = 0;
	std::string fileName = "";

	std::string line;
	std::string pathModel;

	for (const auto& pair : modelTrans) {

		const auto& key = pair.first;
		const auto& vect = pair.second;

		glPushMatrix();

		for (const auto& t : vect) {

			if (t.type == "translate") {
				glTranslatef(t.x, t.y, t.z);
			} else if (t.type == "rotate") {
				glRotatef(t.angle, t.x, t.y, t.z);
			} else if (t.type == "scale") {
				glScalef(t.x, t.y, t.z);
			} else if (t.type == "color") {
				glColor3f(t.x, t.y, t.z);
			}
		}

		fileName = key;
		pathModel = modelsPath + fileName;
		ifstream selectedFile(pathModel);
		float x, y, z;

		if (selectedFile.is_open()) {

			glBegin(GL_TRIANGLES);
			while (std::getline(selectedFile, line)) {

				istringstream iss(line);

				if (iss >> x >> y >> z) {
					glVertex3f(x, y, z);
				}
			}

			glEnd();
			selectedFile.close();
		}
		glPopMatrix();

	}
}

void renderScene(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();
	gluLookAt(posX, posY, posZ,
			  lookAtX, lookAtY, lookAtZ,
			  upX, upY, upZ);

	glPolygonMode(glFaces, glMode);
	
	// Ativar e Desativar os eixos
	if (flagAxis) 
		drawAxis();

	drawModels();

	// End of frame
	glutSwapBuffers();
}

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if (h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load Identity Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set perspective
	gluPerspective(fov, ratio, near, far);

	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

void fillMap(vector<transformation> trans, XMLElement* group) {

	XMLElement* ptrTransform = group->FirstChildElement("transform");

	if (ptrTransform) {

		XMLElement* ptrTRS = ptrTransform->FirstChildElement();

		string type;
		float angle = 0.0f;

		float transX;
		float transY;
		float transZ;

		while (ptrTRS) {

			type = ptrTRS->Name();

			if (type == "rotate") ptrTRS->QueryFloatAttribute("angle", &angle);

			ptrTRS->QueryFloatAttribute("x", &transX);
			ptrTRS->QueryFloatAttribute("y", &transY);
			ptrTRS->QueryFloatAttribute("z", &transZ);

			transformation t(type, angle, transX, transY, transZ);
			trans.push_back(t);

			ptrTRS = ptrTRS->NextSiblingElement();

		}
	}

	XMLElement* ptrNXTGroup = group->FirstChildElement("group");

	while (ptrNXTGroup != nullptr) {
		fillMap(trans, ptrNXTGroup);
		ptrNXTGroup = ptrNXTGroup->NextSiblingElement("group");
	}

	XMLElement* ptrModels = group->FirstChildElement("models");

	if (ptrModels) {

		XMLElement* ptrModel = ptrModels->FirstChildElement("model");
		string modelName = "";

		while (ptrModel) {

			modelName = ptrModel->Attribute("file");
			modelTrans.insert({ modelName, trans });
			ptrModel = ptrModel->NextSiblingElement("model");
		}
	}
}

void loadXML(void) {

	// Criamos um objeto do tipo XMLDocument
	XMLDocument doc;

	if (doc.LoadFile(filePath.c_str()) != tinyxml2::XML_SUCCESS) {
		std::cerr << ">> Error: Reading XML File" << std::endl;
		exit(EXIT_FAILURE);
	}

	XMLElement* ptrRoot = doc.FirstChildElement("world");

	if (!ptrRoot) {
		std::cerr << ">> Error: Node <world> not found!" << std::endl;
		exit(EXIT_FAILURE);
	}

	XMLElement* ptrWindow = ptrRoot->FirstChildElement("window");
	if (ptrWindow) {

		ptrWindow->QueryIntAttribute("width", &width);
		ptrWindow->QueryIntAttribute("height", &height);
	}
	else {
		std::cerr << ">> Error: Node <window> not found" << std::endl;
		exit(EXIT_FAILURE);
	}

	XMLElement* ptrCamera = ptrRoot->FirstChildElement("camera");
	if (ptrCamera) {

		XMLElement* ptrPosition = ptrCamera->FirstChildElement("position");

		if (ptrPosition) {

			ptrPosition->QueryFloatAttribute("x", &posX);
			ptrPosition->QueryFloatAttribute("y", &posY);
			ptrPosition->QueryFloatAttribute("z", &posZ);
		}
		else {
			std::cerr << ">> Error: Node <position> not found" << std::endl;
			exit(EXIT_FAILURE);
		}

		XMLElement* ptrLookAt = ptrCamera->FirstChildElement("lookAt");

		if (ptrLookAt) {

			ptrLookAt->QueryFloatAttribute("x", &lookAtX);
			ptrLookAt->QueryFloatAttribute("y", &lookAtY);
			ptrLookAt->QueryFloatAttribute("z", &lookAtZ);
		}
		else {
			std::cerr << ">> Error: Node <lookAt> not found" << std::endl;
			exit(EXIT_FAILURE);
		}

		XMLElement* ptrUp = ptrCamera->FirstChildElement("up");

		if (ptrUp) {

			ptrUp->QueryFloatAttribute("x", &upX);
			ptrUp->QueryFloatAttribute("y", &upY);
			ptrUp->QueryFloatAttribute("z", &upZ);
		}
		else {
			std::cerr << ">> Error: Node <up> not found" << std::endl;
			exit(EXIT_FAILURE);
		}

		XMLElement* ptrProjection = ptrCamera->FirstChildElement("projection");

		if (ptrProjection) {

			ptrProjection->QueryFloatAttribute("fov", &fov);
			ptrProjection->QueryFloatAttribute("near", &near);
			ptrProjection->QueryFloatAttribute("far", &far);
		}
		else {
			std::cerr << ">> Error: Node <up> not found" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	else {
		std::cerr << ">> Error: Node <camera> not found" << std::endl;
		exit(EXIT_FAILURE);
	}

	alfa = atan2(posZ - lookAtZ, posX - lookAtX);
	beta = atan2(posY - lookAtY, sqrt(pow(posX - lookAtX, 2) + pow(posZ - lookAtZ, 2)));

	// Flag para que não seja necessário reler os modelos caso estejamos a dar reload, basta nos as configurações de camara e posições iniciais
	if (flagReload) return;
	flagReload = true;

	XMLElement* ptrGroup = ptrRoot->FirstChildElement("group");
	vector<transformation> trans;

	while (ptrGroup != nullptr) {
		fillMap(trans, ptrGroup);
		ptrGroup = ptrRoot->NextSiblingElement("group");
	}

	for (const auto& pair : modelTrans) {  // Acessando o par chave-valor
		const auto& key = pair.first;  // Chave (categoria)
		const auto& vec = pair.second;  // Vetor de transformações

		std::cout << "Category: " << key << "\n"; // Imprime a chave
		for (const auto& t : vec) { // Iterando sobre o vetor de transformações
			std::cout << "  Type: " << t.type << ", Angle: " << t.angle
				<< ", Coordinates: (" << t.x << ", " << t.y << ", " << t.z << ")\n";
		}
	}
}

void readKeys(unsigned char c, int x, int y) {

	float radius = sqrt(pow((posX - lookAtX), 2.0f) + pow((posY - lookAtY), 2.0f) + pow((posZ - lookAtZ), 2.0));

	switch (c) {

		case 'r':
			loadXML();
			break;

		case 'd':
			alfa -= 0.1f;
			ExplorerModeCam(radius);
			break;
		case 'a':
			alfa += 0.1f;
			ExplorerModeCam(radius);
			break;
		case 's':
			beta -= 0.1f;
			if (beta < -1.5f)
				beta = -1.5f;
			ExplorerModeCam(radius);
			break;
		case 'w':
			beta += 0.1f;
			if (beta > 1.5f)
				beta = 1.5f;
			ExplorerModeCam(radius);
			break;

		case '+':
				ExplorerModeCam(radius - movementSpeed);
			break;
		case '-':
				ExplorerModeCam(radius + movementSpeed);
			break;

		case '2':
			flagAxis = !flagAxis;
			break;
		case '3':
			if (glFaces == GL_FRONT) {
				glFaces = GL_BACK;
				std::cout << ">> Switched to GL_BACK" << endl;
			}
			else if (glFaces == GL_BACK) {
				glFaces = GL_FRONT_AND_BACK;
				std::cout << ">> Switched to GL_FRONT_AND_BACK" << endl;
			} else {
				glFaces = GL_FRONT;
				std::cout << ">> Switched to GL_FRONT" << endl;
			}
			break;
		case '4':
			if (glMode == GL_LINE) {
				glMode = GL_FILL;
			} else if (glMode == GL_FILL) {
				glMode = GL_POINT;
			} else {
				glMode = GL_LINE;
			}
			break;
	}
	glutPostRedisplay();
}

int main(int argc, char** argv) {

	loadXML();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("CG_TP17_P1");

	glutReshapeFunc(changeSize);
	glutDisplayFunc(renderScene);

	glutKeyboardFunc(readKeys);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	std::cout << "------------------------- Instruções ---------------------------------" << endl;
	std::cout << "											                            " << endl;
	std::cout << "						Cam: ExplorerMode							    " << endl;
	std::cout << "       w: move up / s: move down / a: move left / d: move rigth       " << endl;
	std::cout << "				     +: Zoom In / -: Zoom Out                           " << endl;
	std::cout << "											                            " << endl;
	std::cout << "  r: Reset								                            " << endl;
	std::cout << "						                                                " << endl;
	std::cout << "	2: Turn ON/OFF Axis			                                        " << endl;
	std::cout << "	3: Switch Polygon Faces (GL_FRONT / GL_BACK / GL_FRONT_AND_BACK)    " << endl;
	std::cout << "	4: Switch Polygon Mode (GL_LINE / GL_POINT / GL_FILL)		        " << endl;
	std::cout << "											                            " << endl;
	std::cout << "----------------------------------------------------------------------" << endl;

	glutMainLoop();

	return 1;
}