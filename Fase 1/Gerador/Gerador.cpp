#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>

using namespace std;

std::string generatePlane(float dimensao, int divisoes) {
    std::stringstream triangulos;
    float tamanho_divisoes = dimensao / divisoes;
    float dist_centro = dimensao / 2;
    float x1, x2, z1, z2;
    for (int i = 0; i < divisoes; i++) {
        z1 = (tamanho_divisoes * i) - dist_centro;
        z2 = (tamanho_divisoes * (i + 1)) - dist_centro;

        for (int j = 0; j < divisoes; j++) {
            x1 = (tamanho_divisoes * j) - dist_centro;
            x2 = (tamanho_divisoes * (j + 1)) - dist_centro;

            // Triangulo 1
            triangulos << x1 << ' ' << 0 << ' ' << z1 << '\n';
            triangulos << x1 << ' ' << 0 << ' ' << z2 << '\n';
            triangulos << x2 << ' ' << 0 << ' ' << z2 << '\n';

            // Triangulo 2
            triangulos << x1 << ' ' << 0 << ' ' << z1 << '\n';
            triangulos << x2 << ' ' << 0 << ' ' << z2 << '\n';
            triangulos << x2 << ' ' << 0 << ' ' << z1 << '\n';
        }
    }
    return triangulos.str();
}

std::string generateBox(float dimensao, int divisoes) {
    std::stringstream triangulos;
    float altura = dimensao / 2;
    float tamanho_divisoes = dimensao / divisoes;
    float dist_centro = dimensao / 2;
    float x1, x2, y1, y2, z1, z2;

    //Base e topo
    for (int i = 0; i < divisoes; i++) {
        z1 = (tamanho_divisoes * i) - dist_centro;
        z2 = (tamanho_divisoes * (i + 1)) - dist_centro;

        for (int j = 0; j < divisoes; j++) {
            x1 = (tamanho_divisoes * j) - dist_centro;
            x2 = (tamanho_divisoes * (j + 1)) - dist_centro;

            // y positivo
            // Triangulo 1
            triangulos << x1 << ' ' << altura << ' ' << z1 << '\n';
            triangulos << x1 << ' ' << altura << ' ' << z2 << '\n';
            triangulos << x2 << ' ' << altura << ' ' << z2 << '\n';

            // Triangulo 2
            triangulos << x1 << ' ' << altura << ' ' << z1 << '\n';
            triangulos << x2 << ' ' << altura << ' ' << z2 << '\n';
            triangulos << x2 << ' ' << altura << ' ' << z1 << '\n';

            // y negativo
            // Triangulo 1
            triangulos << x1 << ' ' << -altura << ' ' << z1 << '\n';
            triangulos << x2 << ' ' << -altura << ' ' << z2 << '\n';
            triangulos << x1 << ' ' << -altura << ' ' << z2 << '\n';

            // Triangulo 2
            triangulos << x1 << ' ' << -altura << ' ' << z1 << '\n';
            triangulos << x2 << ' ' << -altura << ' ' << z1 << '\n';
            triangulos << x2 << ' ' << -altura << ' ' << z2 << '\n';
        }
    }

    //Laterais sobre o eixo dos x
    for (int i = 0; i < divisoes; i++) {
        y1 = (tamanho_divisoes * i) - dist_centro;
        y2 = (tamanho_divisoes * (i + 1)) - dist_centro;

        for (int j = 0; j < divisoes; j++) {
            z1 = (tamanho_divisoes * j) - dist_centro;
            z2 = (tamanho_divisoes * (j + 1)) - dist_centro;

            // x negativo
            // Triangulo 1
            triangulos << -dist_centro << ' ' << y1 << ' ' << z1 << '\n';
            triangulos << -dist_centro << ' ' << y1 << ' ' << z2 << '\n';
            triangulos << -dist_centro << ' ' << y2 << ' ' << z2 << '\n';

            // Triangulo 2
            triangulos << -dist_centro << ' ' << y1 << ' ' << z1 << '\n';
            triangulos << -dist_centro << ' ' << y2 << ' ' << z2 << '\n';
            triangulos << -dist_centro << ' ' << y2 << ' ' << z1 << '\n';

            // x positivo
            // Triangulo 1
            triangulos << dist_centro << ' ' << y2 << ' ' << z2 << '\n';
            triangulos << dist_centro << ' ' << y1 << ' ' << z2 << '\n';
            triangulos << dist_centro << ' ' << y1 << ' ' << z1 << '\n';

            // Triangulo 2
            triangulos << dist_centro << ' ' << y2 << ' ' << z2 << '\n';
            triangulos << dist_centro << ' ' << y1 << ' ' << z1 << '\n';
            triangulos << dist_centro << ' ' << y2 << ' ' << z1 << '\n';
        }

    }

    // Laterais sobre o eixo dos z
    for (int i = 0; i < divisoes; i++) {
        y1 = (tamanho_divisoes * i) - dist_centro;
        y2 = (tamanho_divisoes * (i + 1)) - dist_centro;

        for (int j = 0; j < divisoes; j++) {
            x1 = (tamanho_divisoes * j) - dist_centro;
            x2 = (tamanho_divisoes * (j + 1)) - dist_centro;

            // z negativo
            // Triangulo 1
            triangulos << x1 << ' ' << y2 << ' ' << -dist_centro << '\n';
            triangulos << x2 << ' ' << y1 << ' ' << -dist_centro << '\n';
            triangulos << x1 << ' ' << y1 << ' ' << -dist_centro << '\n';

            // Triangulo 2
            triangulos << x1 << ' ' << y2 << ' ' << -dist_centro << '\n';
            triangulos << x2 << ' ' << y2 << ' ' << -dist_centro << '\n';
            triangulos << x2 << ' ' << y1 << ' ' << -dist_centro << '\n';

            // z positivo
            // Triangulo 1
            triangulos << x1 << ' ' << y2 << ' ' << dist_centro << '\n';
            triangulos << x1 << ' ' << y1 << ' ' << dist_centro << '\n';
            triangulos << x2 << ' ' << y1 << ' ' << dist_centro << '\n';

            // Triangulo 2
            triangulos << x1 << ' ' << y2 << ' ' << dist_centro << '\n';
            triangulos << x2 << ' ' << y1 << ' ' << dist_centro << '\n';
            triangulos << x2 << ' ' << y2 << ' ' << dist_centro << '\n';
        }
    }

    return triangulos.str();
}

std::string generateCone(float radius, float height, int slices, int stacks) {
    std::stringstream triangles;
    float angle = 2 * M_PI / slices;
    float heightBase = height / stacks;

    // Base
    for (int i = 0; i < slices; i++) {
        float x1 = radius * sin(angle * i);
        float z1 = radius * cos(angle * i);
        float x2 = radius * sin(angle * (i + 1));
        float z2 = radius * cos(angle * (i + 1));

        triangles << 0 << ' ' << 0 << ' ' << 0 << '\n';
        triangles << x2 << ' ' << 0 << ' ' << z2 << '\n';
        triangles << x1 << ' ' << 0 << ' ' << z1 << '\n';
    }

    for (int i = 0; i < slices; i++) {
        for (int j = 0; j < stacks; j++) {
            float currRadius = radius * (1.0f - float(j) / stacks);
            float nextRadius = radius * (1.0f - float(j + 1) / stacks);

            float x1 = currRadius * sin(angle * i);
            float z1 = currRadius * cos(angle * i);
            float x2 = currRadius * sin(angle * (i + 1));
            float z2 = currRadius * cos(angle * (i + 1));

            float next_x1 = nextRadius * sin(angle * i);
            float next_z1 = nextRadius * cos(angle * i);
            float next_x2 = nextRadius * sin(angle * (i + 1));
            float next_z2 = nextRadius * cos(angle * (i + 1));

            float y = heightBase * j;
            float next_y = heightBase * (j + 1);


            // Left Triangles
            triangles << x1 << ' ' << y << ' ' << z1 << '\n';
            triangles << x2 << ' ' << y << ' ' << z2 << '\n';
            triangles << next_x1 << ' ' << next_y << ' ' << next_z1 << '\n';

            // Right Triangles
            triangles << x2 << ' ' << y << ' ' << z2 << '\n';
            triangles << next_x2 << ' ' << next_y << ' ' << next_z2 << '\n';
            triangles << next_x1 << ' ' << next_y << ' ' << next_z1 << '\n';

        }
    }

    return triangles.str();
}


std::string generateSphere(float radius, int slices, int stacks) {
    std::stringstream triangles;
    float dPhi = 2 * M_PI / slices;
    float dTheta = M_PI / stacks;

    for (int j = 0; j < stacks; j++) {
        float angle = j * dTheta;
        float next_angle = (j + 1) * dTheta;

        for (int i = 0; i < slices; i++) {
            float space = dPhi * i;
            float next_space = dPhi * (i + 1);

            float x1 = radius * sin(angle) * cos(space);
            float z1 = radius * sin(angle) * sin(space);

            float x2 = radius * sin(angle) * cos(next_space);
            float z2 = radius * sin(angle) * sin(next_space);

            float y = radius * cos(angle);

            float next_x1 = radius * sin(next_angle) * cos(space);
            float next_z1 = radius * sin(next_angle) * sin(space);

            float next_x2 = radius * sin(next_angle) * cos(next_space);
            float next_z2 = radius * sin(next_angle) * sin(next_space);

            float next_y = radius * cos(next_angle);

            // Left Triangles
            triangles << x1 << ' ' << y << ' ' << z1 << '\n';
            triangles << x2 << ' ' << y << ' ' << z2 << '\n';
            triangles << next_x1 << ' ' << next_y << ' ' << next_z1 << '\n';

            // Right Triangles
            triangles << x2 << ' ' << y << ' ' << z2 << '\n';
            triangles << next_x2 << ' ' << next_y << ' ' << next_z2 << '\n';
            triangles << next_x1 << ' ' << next_y << ' ' << next_z1 << '\n';
        }
    }
    return triangles.str();
}

std::string generateCylinder(float radius, float height, int slices) {
    std::stringstream triangles;
    float angleStep = 2.0f * M_PI / slices;
    for (int i = 0; i < slices; i++) {
        float angle = i * angleStep;
        float next_angle = (i + 1) * angleStep;
        float x1 = radius * sin(angle);
        float z1 = radius * cos(angle);
        float x2 = radius * sin(next_angle);
        float z2 = radius * cos(next_angle);
        
        // Bottom Base
        triangles << x2 << ' ' << 0 << ' ' << z2 << '\n';
        triangles << x1 << ' ' << 0 << ' ' << z1 << '\n';
        triangles << 0 << ' ' << 0 << ' ' << 0 << '\n';

        // Top Base
        triangles << x1 << ' ' << height << ' ' << z1 << '\n';
        triangles << x2 << ' ' << height << ' ' << z2 << '\n';
        triangles << 0 << ' ' << height << ' ' << 0 << '\n';

        // Left Triangle
        triangles << x1 << ' ' << 0 << ' ' << z1 << '\n';
        triangles << x2 << ' ' << 0 << ' ' << z2 << '\n';
        triangles << x1 << ' ' << height << ' ' << z1 << '\n';

        // Right Triangle
        triangles << x2 << ' ' << 0 << ' ' << z2 << '\n';
        triangles << x2 << ' ' << height << ' ' << z2 << '\n';
        triangles << x1 << ' ' << height << ' ' << z1 << '\n';
    }
    return triangles.str();
}

std::string generateTorus(float outerRadious, float innerRadius, int rings, int stacks){
    std::stringstream triangles;
    float outerAngle = 2 * M_PI / rings;
    float innerAngle = 2 * M_PI / stacks;
    float x1, x2, nextx1, nextx2, y1, y2, z1, z2, nextz1, nextz2;
    
    for(int i = 0; i < rings; i++){
        for(int j = 0; j < stacks; j++){
            // x e z do primeiro ponto do atual "anel"
            x1 = (outerRadious + innerRadius * cos(innerAngle * j)) * sin(outerAngle * i);
            z1 = (outerRadious + innerRadius * cos(innerAngle * j)) * cos(outerAngle * i);
            
            // x e z do segundo ponto do atual "anel"
            x2 = (outerRadious + innerRadius * cos(innerAngle * (j + 1))) * sin(outerAngle * i);
            z2 = (outerRadious + innerRadius * cos(innerAngle * (j + 1))) * cos(outerAngle * i);
            
            // x e z do primeiro ponto do proximo "anel"
            nextx1 = (outerRadious + innerRadius * cos(innerAngle * j)) * sin(outerAngle * (i + 1));
            nextz1 = (outerRadious + innerRadius * cos(innerAngle * j)) * cos(outerAngle * (i + 1));
            
            // x e z do segundo ponto do proximo "anel"
            nextx2 = (outerRadious + innerRadius * cos(innerAngle * (j + 1))) * sin(outerAngle * (i + 1));
            nextz2 = (outerRadious + innerRadius * cos(innerAngle * (j + 1))) * cos(outerAngle * (i + 1));
            
            // y dos priemiros pontos de cada anel
            y1 = innerRadius * sin(innerAngle * j);
            // y dos segundos pontos de cada anel
            y2 = innerRadius * sin(innerAngle * (j + 1));

            // Triangulo 1
            triangles << x1 << ' ' << y1 << ' ' << z1 << '\n';
            triangles << nextx2 << ' ' << y2 << ' ' << nextz2 << '\n';
            triangles << x2 << ' ' << y2 << ' ' << z2 << '\n';

            // Triangulo 2
            triangles << x1 << ' ' << y1 << ' ' << z1 << '\n';
            triangles << nextx1 << ' ' << y1 << ' ' << nextz1 << '\n';
            triangles << nextx2 << ' ' << y2 << ' ' << nextz2 << '\n';
        }
    }
    return triangles.str();
}

std::string gera_ficheiro(int argc, char const* argv[]) {
    std::string ficheiro = "";
    if (std::string(argv[1]) == "sphere") {
        ficheiro = generateSphere(std::stof(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]));
    }
    else if (std::string(argv[1]) == "box") {
        ficheiro = generateBox(std::stof(argv[2]), std::stoi(argv[3]));
    }
    else if (std::string(argv[1]) == "cone") {
        ficheiro = generateCone(std::stof(argv[2]), std::stof(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
    }
    else if (std::string(argv[1]) == "plane") {
        ficheiro = generatePlane(std::stof(argv[2]), std::stoi(argv[3]));
    }
    else if (std::string(argv[1]) == "cylinder") {
        ficheiro = generateCylinder(std::stof(argv[2]), std::stof(argv[3]), std::stoi(argv[4]));
    }
    else if (std::string(argv[1]) == "torus") {
        ficheiro = generateTorus(std::stof(argv[2]), std::stof(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
    }
    return ficheiro;
}



int isFloat(char const* n) {
    int r = 1;
    int p = 0;
    for (int i = 0; i < strlen(n); i++) {
        if (!isdigit(n[i])) {
            if (i != 0 && n[i] == '.' && p == 0) {
                p = 1;
            }
            else {
                r = 0;
            }
        }
    }
    return r;
}


int isInt(char const* n) {
    int r = 1;
    for (int i = 0; i < strlen(n); i++) {
        if (!isdigit(n[i])) {
            r = 0;
        }
    }
    return r;
}


int valida_input(int argc, char const* argv[]) {
    int r = 0;
    if ((std::string(argv[1]) == "sphere" && argc == 6) || (std::string(argv[1]) == "box" && argc == 5) ||
        (std::string(argv[1]) == "cone" && argc == 7) || (std::string(argv[1]) == "plane" && argc == 5) ||
        (std::string(argv[1]) == "cylinder" && argc == 6) || (std::string(argv[1]) == "torus" && argc == 7)) {
        r = 1;

        for (int i = 2; i < argc - 1; i++) {
            if (!isInt(argv[i]) && !isFloat(argv[i])) {
                r = 0;
            }
        }
    }
    return r;

}


int main(int argc, char const* argv[]) {

    if (valida_input(argc, argv)) {
        ofstream ficheiro;
        ficheiro.open("..\\models\\" + std::string(argv[argc - 1]));

        ficheiro << gera_ficheiro(argc, argv);

        ficheiro.close();
    }
    else {
        cout << "Input invalido" << endl;
    }

    return 0;
}