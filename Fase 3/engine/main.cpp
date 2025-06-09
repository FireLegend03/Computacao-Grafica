#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
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

// Advanced camera settings

bool camModeFPS = true;

float lastX = -1.0f;
float lastY = -1.0f;
float startX = 0.0f;
float startY = 0.0f;
float radius = 0.0f;

float camYaw = 0.0f;
float camPitch = 0.0f;

bool dragging = false;
int tracking = 0;

// Paths
std::string filePath = "..\\..\\xml\\sistema_solar.xml";
std::string modelsPath = "..\\..\\models\\";

// Configs
bool mousePressed = false;
float movementSpeed = 1;

bool flagCatMull = true;
bool flagAxis = true;
bool flagReload = false;
int glMode = GL_LINE;
int glFaces = GL_FRONT; 

// VBO's
GLuint vertices;

struct transformation {

	std::string type;
	float angle;
	float total_time;
	std::string align;
	std::vector<float> points;
	float x;
	float y;
	float z;

	transformation (std::string type, float angle, float total_time, std::string align,
				   std::vector<float> points, float x, float y, float z)
					: type(type), angle(angle), total_time(total_time), align(align),
					 points(points), x(x), y(y), z(z) {
	}
};

std::multimap<string, std::vector<transformation>> modelTrans;
std::vector<int> nPoints;

// time Control
bool flagBackTrack = false;
bool flagPause = false;
double current_time = 0.f;
float pausedInterval = 0.f;
float initBackTrack = 0.f;
float initPause = 0.f;
float timeLost = 0.f;
double lastFrameTime = 0.f;


void ExplorerModeCam(float radius) {

	posX = radius * cos(beta) * cos(alfa);
	posY = radius * sin(beta);
	posZ = radius * cos(beta) * sin(alfa);
}
//
//void FPSModeCam() {
//
//	float radYaw = - camYaw * M_PI / 180.0f;
//	float radPitch = camPitch * M_PI / 180.0f;
//
//	lookAtX = posX + cos(radPitch) * sin(radYaw);
//	lookAtY = posY + sin(radPitch);
//	lookAtZ = posZ - cos(radPitch) * cos(radYaw);
//}

void drawAxis(void) {

	glBegin(GL_LINES);

	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(800.0f, 0.0f, 0.0f);

	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 800.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 800.0f);

	glColor3f(1.0f, 1.0f, 1.0f);

	glEnd();
}

// Fun��o auxiliar para multiplicar uma matriz por um vetor
void multMatrixVector(float matrix[4][4], float* vector, float* result) {
	for (int j = 0; j < 4; ++j)
	{
		result[j] = 0;
		for (int k = 0; k < 4; ++k)
		{
			result[j] += vector[k] * matrix[j][k];
		}
	}
}

void getCatmullRomPoint(float t, float* p0, float* p1, float* p2, float* p3, float* pos, float* deriv) {

	// catmull-rom matrix
	float M[4][4] = { {-0.5f,  1.5f, -1.5f,  0.5f},
					  { 1.0f, -2.5f,  2.0f, -0.5f},
					  {-0.5f,  0.0f,  0.5f,  0.0f},
					  { 0.0f,  1.0f,  0.0f,  0.0f} };

	for (int i = 0; i < 3; i++) {
		float P[4] = { p0[i], p1[i], p2[i], p3[i] };
		float A[4];
		
		//  M * P = A
		multMatrixVector(M, P, A);

		// calculamos a posição, pos = T * A
		pos[i] = powf(t, 3) * A[0] + powf(t, 2) * A[1] + t * A[2] + A[3];

		// calculamos a derivada,  deriv = T' * A
		deriv[i] = 3 * powf(t, 2) * A[0] + 2 * t * A[1] + A[2];
	}

}

// returns the point in the curve
void getGlobalCatmullRomPoint(float time, float* pos, float* deriv, std::vector<float> points) {

	int numPoints = points.size() / 3; // numero de pontos na curva
	// Calcular o segmento atual da nossa curva
	float t = time * ((float)points.size() / 3); 
	int index = floor(t);  //indice do segmento atual
	t = t - index; // tempo dentro do segmento atual

	// indices store the points
	int indices[4];
	indices[0] = (index + numPoints - 1) % numPoints;
	indices[1] = (indices[0] + 1) % numPoints;
	indices[2] = (indices[1] + 1) % numPoints;
	indices[3] = (indices[2] + 1) % numPoints;

	// Matriz que armazena as coordenadas dos 4 pontos atuais
	float p[4][3];
	for (int i = 0; i < 4; i++) {
		p[i][0] = points[indices[i] * 3];
		p[i][1] = points[indices[i] * 3 + 1];
		p[i][2] = points[indices[i] * 3 + 2];
	}

	getCatmullRomPoint(t, p[0], p[1], p[2], p[3], pos, deriv);
}

// Função que desenha a curva de catmull-rom
void catmullRomCurve(std::vector<float> points) {
	float pos[3];
	float deriv[3];

	glBegin(GL_LINE_LOOP);
	glColor3f(1.0f, 1.0f, 1.0f);

	for (float t = 0; t < 1; t += 0.01) {
		getGlobalCatmullRomPoint(t, pos, deriv, points);
		glVertex3f(pos[0], pos[1], pos[2]);
	}
	glEnd();
}

void buildRotMatrix(float* x, float* y, float* z, float* m) {

	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}

void cross(float* a, float* b, float* res) {

	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

void normalize(float* a) {

	float l = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0] / l;
	a[1] = a[1] / l;
	a[2] = a[2] / l;
}

void alignCatMull(float *deriv, float *m) {

	float x[3] = { deriv[0],deriv[1],deriv[2] };
	float y[3];
	float z[3];
	float upVec[3] = { 0,1,0 };

	normalize(x);
	cross(x, upVec, z);
	normalize(z);
	cross(z, x, y);
	normalize(y);
	buildRotMatrix(x, y, z, m);
	glMultMatrixf(m);
	
}

void catmullRom(transformation& t) {
	float pos[3];
	float deriv[3];
	
	if (flagCatMull) catmullRomCurve(t.points);

	// dividir o tempo atual pelo tempo total para obter o tempo entre os valores 0 e 1
	// para depois utilizar no calculo do segmento atual da nossa curva
	float time = current_time / t.total_time; 

	getGlobalCatmullRomPoint(time, pos, deriv, t.points);

	// Atualizar a posição do modelo
	glTranslatef(pos[0], pos[1], pos[2]);
	if (t.align == "true") {
		float m[16];
		alignCatMull(deriv, m);
	}
}

void getTimedAngle(transformation& t) {
	float rotation_speed = 360 / t.total_time;
	t.angle = current_time * rotation_speed;
}

void drawModels(void) {

	int modelIndex = 0;
	int initVertex = 0;
	std::string fileName = "";

	std::string pathModel;

	glBindBuffer(GL_ARRAY_BUFFER, vertices);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	for (auto& pair : modelTrans) {

		auto& key = pair.first;
		auto& transformations = pair.second;

		glPushMatrix();

		for (auto& t : transformations) {

			if (t.type == "translate") {
				glTranslatef(t.x, t.y, t.z);
			} else if (t.type == "rotate") {
				glRotatef(t.angle, t.x, t.y, t.z);
			}
			else if (t.type == "timedRotate") {
				getTimedAngle(t);
				glRotatef(t.angle, t.x, t.y, t.z);
			} else if (t.type == "scale") {
				glScalef(t.x, t.y, t.z);
			} else if (t.type == "color") {
				glColor3f(t.x, t.y, t.z);
			} else if (t.type == "catmull") {
				catmullRom(t);
			}
		}

		fileName = key;
		pathModel = modelsPath + fileName;
		ifstream selectedFile(pathModel);
		float x, y, z;

		int pointsPerModel = nPoints[modelIndex];
		glDrawArrays(GL_TRIANGLES, initVertex, pointsPerModel);
		initVertex += pointsPerModel;
		modelIndex++;
		
		glPopMatrix();
	}
}

void updateTime() {
	if (!flagPause) {
		double now = glutGet(GLUT_ELAPSED_TIME) / 1000.0 - (timeLost / 1000.0);
		double deltaTime = now - lastFrameTime;
		lastFrameTime = now;

		if (!flagBackTrack) {
			current_time += deltaTime;
		}
		else {
			current_time -= deltaTime;
			if (current_time < 0.0) current_time = 0.0;
		}
	}
}

void renderScene(void) {

	updateTime();

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

void prepareData() {
	// Criar um vetor com as coordenadas dos pontos
	vector<float> points;

	GLuint pointCounter = 0;
	std::string fileName = "";

	std::string line;
	std::string pathModel;

	for ( auto& pair : modelTrans) {

		pointCounter = 0;
		auto& key = pair.first;
		auto& vect = pair.second;

		fileName = key;
		pathModel = modelsPath + fileName;
		ifstream selectedFile(pathModel);
		float x, y, z;

		if (selectedFile.is_open()) {

			while (std::getline(selectedFile, line)) {

				istringstream iss(line);

				if (iss >> x >> y >> z) {
					points.push_back(x);
					points.push_back(y);
					points.push_back(z);
					pointCounter++;
				}
			}

			selectedFile.close();
		}

		nPoints.push_back(pointCounter);
	}

	// Criar o VBO
	glEnableClientState(GL_VERTEX_ARRAY);
	glGenBuffers(1, &vertices);

	// Copiar o vetor para a memória gráfica
	glBindBuffer(GL_ARRAY_BUFFER, vertices);
	glBufferData(
		GL_ARRAY_BUFFER, 
		sizeof(float) * points.size(), 
		points.data(),
		GL_STATIC_DRAW); 
}

void fillMap(vector<transformation> trans, XMLElement* group) {

	XMLElement* ptrTransform = group->FirstChildElement("transform");

	if (ptrTransform) {

		XMLElement* ptrTRS = ptrTransform->FirstChildElement();

		string type;
		float angle = 0.0f;
		float total_time = 0.0f;
		float current_time = 0.0f;
		const char* align = "false";
		vector<float> points;
		float transX = 0.0f;
		float transY = 0.0f;
		float transZ = 0.0f;

		while (ptrTRS) {

			type = ptrTRS->Name();
			if (type == "translate") {
				XMLElement* ptrPoints = ptrTRS->FirstChildElement("point");
				if (ptrPoints) { //Translate de catmull
					while (ptrPoints) {
						//ptrTRS->QueryFloatAttribute("align", &align); tratar do align mais tarde
						ptrPoints->QueryFloatAttribute("x", &transX);
						ptrPoints->QueryFloatAttribute("y", &transY);
						ptrPoints->QueryFloatAttribute("z", &transZ);

						points.push_back(transX);
						points.push_back(transY);
						points.push_back(transZ);

						ptrPoints = ptrPoints->NextSiblingElement("point");
					}
					type = "catmull";
					ptrTRS->QueryFloatAttribute("time", &total_time);
					ptrTRS->QueryStringAttribute("align", &align);

				} else { //Translate normal
					ptrTRS->QueryFloatAttribute("x", &transX);
					ptrTRS->QueryFloatAttribute("y", &transY);
					ptrTRS->QueryFloatAttribute("z", &transZ);
				}
			} else if (type == "rotate") { // Rotate normal
				if (ptrTRS->Attribute("time") != nullptr) {
					type = "timedRotate";
					ptrTRS->QueryFloatAttribute("time", &total_time);
				}
				ptrTRS->QueryFloatAttribute("angle", &angle);
				ptrTRS->QueryFloatAttribute("x", &transX);
				ptrTRS->QueryFloatAttribute("y", &transY);
				ptrTRS->QueryFloatAttribute("z", &transZ);
			}
			else {
				ptrTRS->QueryFloatAttribute("x", &transX);
				ptrTRS->QueryFloatAttribute("y", &transY);
				ptrTRS->QueryFloatAttribute("z", &transZ);
			}

			transformation t(type, angle, total_time, align, points, transX, transY, transZ);
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
	float dx = lookAtX - posX;
	float dy = lookAtY - posY;
	float dz = lookAtZ - posZ;

	radius = sqrtf(dx * dx + dy * dy + dz * dz);

	// Flag para que não seja necessário reler os modelos caso estejamos a dar reload, basta nos as configurações de camara e posições iniciais
	if (flagReload) return;
	flagReload = true;

	XMLElement* ptrGroup = ptrRoot->FirstChildElement("group");
	vector<transformation> trans;

	while (ptrGroup != nullptr) {
		fillMap(trans, ptrGroup);
		ptrGroup = ptrRoot->NextSiblingElement("group");
	}
}

void releaseKeys(unsigned char c, int x, int y) {

	switch (c) {
		case 'b':
			if (flagBackTrack && !flagPause) {
				flagBackTrack = false;
				timeLost += glutGet(GLUT_ELAPSED_TIME) - initBackTrack;
				lastFrameTime = glutGet(GLUT_ELAPSED_TIME) / 1000.0 - timeLost / 1000.0;
			}
			break;
	}
}

void readKeys(unsigned char c, int x, int y) {

	float radius = sqrt(pow((posX - lookAtX), 2.0f) + pow((posY - lookAtY), 2.0f) + pow((posZ - lookAtZ), 2.0));

	switch (c) {

		// option keys
	case 'r':
		loadXML();
		current_time = 0.f;
		break;
	case 'p':
		if (!flagBackTrack) {
			if (!flagPause) {
				flagPause = true;
				initPause = glutGet(GLUT_ELAPSED_TIME);
			}
			else {
				timeLost += glutGet(GLUT_ELAPSED_TIME) - initPause;
				flagPause = false;
			}
		}
		break;
	case 'b':
		if (!flagPause)
			if (!flagBackTrack) {
				flagBackTrack = true;
				initBackTrack = glutGet(GLUT_ELAPSED_TIME);
				lastFrameTime = glutGet(GLUT_ELAPSED_TIME) / 1000.0 - timeLost / 1000.0;
			}
		break;

		// TODO: se implementar o movimento da camara só com o rato, remover estas antigas opções
		
		// movement keys
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


		// camera keys
	case '+':
		ExplorerModeCam(radius - movementSpeed);
		break;
	case '-':
		ExplorerModeCam(radius + movementSpeed);
		break;

		// Draw Option keys
	case '1':
		if (camModeFPS) {
			lookAtX = lookAtY = lookAtZ = 0.f;
		}


		camModeFPS = !camModeFPS;
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
		}
		else {
			glFaces = GL_FRONT;
			std::cout << ">> Switched to GL_FRONT" << endl;
		}
		break;
	case '4':
		if (glMode == GL_LINE) {
			glMode = GL_FILL;
		}
		else if (glMode == GL_FILL) {
			glMode = GL_POINT;
		}
		else {
			glMode = GL_LINE;
		}
		break;
	case '5':
		flagCatMull = !flagCatMull;
		break;

	}
	glutPostRedisplay();
}

/*
void readMouseButton(int button, int state, int xx, int yy) {
	
	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
		else
			tracking = 0;
	}
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alfa += (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2) {

			radius -= yy - startY;
			if (radius < 3)
				radius = 3.0;
		}
		tracking = 0;
	}
}

void readMouseMotion(int xx, int yy) {
	int deltaX, deltaY;
	int alphaAux, betaAux;
	int rAux;

	if (!tracking)
		return;

	deltaX = xx - startX;
	deltaY = yy - startY;

	if (tracking == 1) {

		alphaAux = alfa - deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0)
			betaAux = 85.0;
		else if (betaAux < -85.0)
			betaAux = -85.0;

		rAux = radius;
	} else if (tracking == 2) {

		alphaAux = alfa;
		betaAux = beta;
		rAux = radius - deltaY;
		if (rAux < 3)
			rAux = 3;
	}
	posX = rAux * sin(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
	posZ = rAux * cos(alphaAux * 3.14 / 180.0) * cos(betaAux * 3.14 / 180.0);
	posY = rAux * sin(betaAux * 3.14 / 180.0);
}
*/


void printInfo() {
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
	std::cout << "	5: Turn ON/OFF CatMull-Rom Curves                                   " << endl;
	std::cout << "											                            " << endl;
	std::cout << "----------------------------------------------------------------------" << endl;
}

int main(int argc, char** argv) {

	lastFrameTime = glutGet(GLUT_ELAPSED_TIME) / 1000.0;

	loadXML();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow("CG_TP17_P1");

	glutReshapeFunc(changeSize);
	glutDisplayFunc(renderScene);
	glutIdleFunc(renderScene);

	glutKeyboardFunc(readKeys);
	glutKeyboardUpFunc(releaseKeys);
	//glutMouseFunc(readMouseButton);
	//glutMotionFunc(readMouseMotion);

	GLenum err = glewInit();
	if (GLEW_OK != err) {
		fprintf(stderr, ">> Error: %s\n", glewGetErrorString(err));
		return 1;
	}
	fprintf(stdout, ">> Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	printInfo();
	prepareData();

	glutMainLoop();

	return 1;
}