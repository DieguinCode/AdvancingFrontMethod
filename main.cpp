#include <vector>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "alg.hpp"
#include "objalgorithms.hpp"
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Function to create and load a shader
static GLuint createShaderProgram(const char* vertexShaderSource, const char* fragmentShaderSource) {
    // Compiling vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    GLint success;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "Erro ao compilar o vertex shader:\n" << infoLog << std::endl;
    }

    // Compiling fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "Erro ao compilar o fragment shader:\n" << infoLog << std::endl;
    }

    // Link shader
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "Erro ao linkar o programa de shader:\n" << infoLog << std::endl;
    }

    // Delete mid shaders
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}

// Function to draw input's points
static void drawInputPoints(const std::vector<vec2>& points, GLuint shaderProgram) {
    // Vertex Array Object (VAO)
    GLuint VAO;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    // Vertex Buffer Object (VBO)
    GLuint VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(vec2), points.data(), GL_STATIC_DRAW);

    // Attributes's layout
    GLint posAttrib = glGetAttribLocation(shaderProgram, "aPos");
    glVertexAttribPointer(posAttrib, 2, GL_DOUBLE, GL_FALSE, sizeof(vec2), (void*)0);
    glEnableVertexAttribArray(posAttrib);

    glUseProgram(shaderProgram);
    glPointSize(5.0f); // Size

    glDrawArrays(GL_POINTS, 0, points.size());

    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}

// Function to draw Convex Hull's Edges
static void drawConvexHull(const std::vector<vec2>& points, GLuint shaderProgram) {
    // Check if there are enough points to build a convex hull
    if (points.size() < 3) {
        std::cerr << "� necess�rio pelo menos 3 pontos para desenhar o convex hull!" << std::endl;
        return;
    }

    // Vertex Array Object (VAO)
    GLuint VAO;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    // Vertex Buffer Object (VBO) 
    GLuint VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(vec2), points.data(), GL_STATIC_DRAW);

    // Attributes's layout
    GLint posAttrib = glGetAttribLocation(shaderProgram, "aPos");
    glVertexAttribPointer(posAttrib, 2, GL_DOUBLE, GL_FALSE, sizeof(vec2), (void*)0);
    glEnableVertexAttribArray(posAttrib);

    glUseProgram(shaderProgram);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Draw just the lines

    glDrawArrays(GL_LINE_LOOP, 0, points.size());

    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}

static void drawTriangulation(const std::vector<vec2>& points, GLuint shaderProgram) {
    
    if (points.size() < 3) {
        std::cerr << "Erro: Menos que 3 pontos para desenhar um triangulo!" << std::endl;
        return;
    }

    // Vertex Array Object (VAO)
    GLuint VAO;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    // Vertex Buffer Object (VBO) 
    GLuint VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(vec2), points.data(), GL_STATIC_DRAW);

    // Attributes's layout
    GLint posAttrib = glGetAttribLocation(shaderProgram, "aPos");
    glVertexAttribPointer(posAttrib, 2, GL_DOUBLE, GL_FALSE, sizeof(vec2), (void*)0);
    glEnableVertexAttribArray(posAttrib);

    glUseProgram(shaderProgram);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Draw just the lines

    glDrawArrays(GL_TRIANGLES, 0, points.size());

    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}

// Vertex Shader Source
const char* vertexShaderSource = R"(
#version 330 core
layout(location = 0) in vec2 aPos;
uniform vec2 minPos;
uniform vec2 maxPos;

void main() {
    vec2 scaledPos = (aPos - minPos) / (maxPos - minPos) * 2.0 - 1.0;
    gl_Position = vec4(scaledPos.x, scaledPos.y, 0.0, 1.0);
}
)";

// Fragment Shader Source
const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;
void main() {
    FragColor = vec4(1.0, 1.0, 1.0, 1.0); // Cor branca
}
)";

int main() {
    if (!glfwInit()) {
        std::cerr << "Falha ao inicializar GLFW!" << std::endl;
        return -1;
    }


    GLFWwindow* window = glfwCreateWindow(800, 600, "Janela OpenGL", nullptr, nullptr);
    if (!window) {
        std::cerr << "Falha ao criar janela GLFW!" << std::endl;
        glfwTerminate();
        return -1;
    }


    glfwMakeContextCurrent(window);


    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Falha ao inicializar GLAD!" << std::endl;
        return -1;
    }


    glViewport(0, 0, 800, 600);

    // Set the BackGround Color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black

    // Here would be possible to read some .obj input file

    //std::vector<vec2> inputPoints = {
        //{1, 1}, {2, 1}, {4, 3}, {3, 2}, {1, -3}, {3, -2}, {-2, -1}, {-4, -3}, {-3, 3}, {-2, 2}
    //};

    //std::vector<vec2> inputPoints = {
    //    {1.1, 1.2}, {2.3, 1.4}, {4.5, 3.6}, {3.7, 4.01}, {1.9, -3.01}, {3.12, -2.34},
    //        {-2.56, -1.78}, {-4.91, -3.123}, {-3.456, 3.789}, {-2.1111, 2.2222}
    //};
    
    //vec2 a{ 1.0, 1.0 };
    //vec2 b{ 1.0, -1.0 };

    //bool tmp_cond = a.is_left(b);
    //if (tmp_cond) { std::cout << "'A' est� � esquerda de 'B'" << std::endl; }
    //else { std::cout << "'A' est� � direita de 'B'" << std::endl; }

    //SQUARE
    //std::vector<vec2> inputPoints{
    //    {1, 2}, {1, 7}, {6, 2}, {6, 7}
    //    , {3.5, 4.5}  //Centroid 
    //};

    //BETTER CIRCLE CENTROID
    //std::vector<vec2> inputPoints{
    //    {0,0}, {0,-1}, {0,1}, {1, 0}, {-1, 0},
    //    {2,2}, {-2,-2}, {2,-2}, {-2,2},
    //    {3,0}, {0,3}, {-3, 0}, {0, -3}
    //};

     random_device rd;
     mt19937 gen(rd());
     uniform_real_distribution<> dis(-500.0, 500.0);
     vector<vec2> inputPoints{};
     int iterationSize = 20; //100 or 1000
     int pointCount = 0;
     for (int i = 0; i < iterationSize; i++) {
         //vec2 point = vec2(round(dis(gen)), round(dis(gen)));
         vec2 point = vec2(dis(gen), dis(gen));
         inputPoints.push_back(point);
         //inputPoints.push_back(vec2(((dis(gen))), ((dis(gen)))));
         pointCount += 1;
         //cout << point << endl;
         if (pointCount >= iterationSize) break;
     }
     //cout << endl;

     vector<vec2> copy{ inputPoints.begin(), inputPoints.end() };

     std::vector<vec2> convexHull = jarvis(&copy);
    // vector<vec2> triangulation = adf2(inputPoints, convexHull);

    vector<vec2> triangulation = vector<vec2>();

    Mesh* onix = ObjUtils::loadMesh("onix.obj");
    vector<vector<vec2>> onixADFResult = advancingFrontObjFirstStep(onix);
    for (vector<vec2> subObjectPoints : onixADFResult) {
        for (vec2 point : subObjectPoints) {
            triangulation.push_back(point);
        }
    }
    Mesh* triangulatedOnix = advancingFrontObjSecondStep(onixADFResult, "triangulated_onix");
    ObjUtils::saveMesh(triangulatedOnix, "triangulated_onix.obj");
    ObjUtils::saveMesh(onix, "onix_nochange.obj");


    // Build the shader program
    GLuint shaderProgram = createShaderProgram(vertexShaderSource, fragmentShaderSource);

    // Calculate the limits (points)
    vec2 minPos = triangulation[0];
    vec2 maxPos = triangulation[0];
    for (const auto& point : triangulation) {
        minPos.x = std::min(minPos.x, point.x);
        minPos.y = std::min(minPos.y, point.y);
        maxPos.x = std::max(maxPos.x, point.x);
        maxPos.y = std::max(maxPos.y, point.y);
    }

    // Add a padding to improve visualization
    double padding = 0.1 * std::max(maxPos.x - minPos.x, maxPos.y - minPos.y);
    minPos.x -= padding;
    minPos.y -= padding;
    maxPos.x += padding;
    maxPos.y += padding;

    while (!glfwWindowShouldClose(window)) {
        
        glClear(GL_COLOR_BUFFER_BIT);

        // To define Shader's uniforms
        glUseProgram(shaderProgram);
        GLint minPosLoc = glGetUniformLocation(shaderProgram, "minPos");
        glUniform2f(minPosLoc, minPos.x, minPos.y);
        GLint maxPosLoc = glGetUniformLocation(shaderProgram, "maxPos");
        glUniform2f(maxPosLoc, maxPos.x, maxPos.y);

        // Draw the points
        drawInputPoints(inputPoints, shaderProgram);

        // Draw the edges
        drawConvexHull(convexHull, shaderProgram);

        //Draw the Triangulation
        drawTriangulation(triangulation, shaderProgram);

        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;
}
