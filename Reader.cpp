#include "Reader.h"
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>

Reader::Reader() {}

Reader::Reader(const Reader& ) {}

Reader::~Reader() {}

/**
 * @brief Reader::readMeshFile reads file with 2D coordinates of vertices and indices linking them
 * @param mesh object in which information is put
 * @param path file name
 * @return
 */
bool Reader::readMeshFile(Mesh<Vertex2D, Triangle> &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++) {
        double x, y;

        input >> x;
        input >> y;
        if (input.eof())
            break;

        Vertex2D v = Vertex2D(x, y);
        mesh.addVertex(v);
    }

    int num_topSimplexes;
    input >> num_topSimplexes;

    if(num_topSimplexes == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topSimplexes);

    //legge i top simplessi
    for (int i = 0; i < num_topSimplexes; i++) {
        int v[3];
        for (int j = 0; j < 3; j++)
            input >> v[j];
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.addTopSimplex(t);
    }

    return true;
}

/**
 * @brief Reader::readMeshFile reads file with 3D coordinates of vertices and indices linking them
 * @param mesh object in which information is put
 * @param path file name
 * @return
 */
bool Reader::readMeshFile(Mesh<Vertex3D, Triangle> &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //reads vertices
    for (int i = 0; i < num_vertices; i++) {
        double x, y, z;

        input >> x;
        input >> y;
        input >> z;
        if (input.eof())
            break;

        Vertex3D v = Vertex3D(x, y, z);
        mesh.addVertex(v);
    }

    int num_topSimplexes;
    input >> num_topSimplexes;

    if(num_topSimplexes == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topSimplexes);

    //reads top simplexes
    for (int i = 0; i < num_topSimplexes; i++) {
        int v[3];
        for (int j = 0; j < 3; j++)
            input >> v[j];
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.addTopSimplex(t);
    }

    return true;
}

/**
 * @author Giulia Picciau
 * @date 04-15-2014
 * @brief Reader::readOFFMesh reads triangle meshes in .off format
 * @param mesh where to put the information
 * @param path path to file
 * @return
 */
bool Reader::readOFFMesh(Mesh<Vertex3D, Triangle> &mesh, string path){

    ifstream input(path.c_str());

    if(input.is_open()==false){
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string line;

    getline(input,line);
    if(line.compare("OFF")){
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }
    getline(input,line);
    istringstream iss(line);
    //iss.
    unsigned long int num_vertices, num_topsimplexes, num_edges; //edges is useless

    iss >> num_vertices;//input >> num_vertices;
    iss >> num_topsimplexes;//input >> num_topsimplexes;
    iss >> num_edges;

    cout<<"V "<<num_vertices<<" T "<<num_topsimplexes<<endl;

    if(num_vertices == 0 || num_topsimplexes == 0){
        cerr<< "Number of simplexes is 0 Not a valid .off file "<<path<<endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //insert vertices
    for(unsigned long int i=0;i<num_vertices;i++){
        double x,y,z;
        input >> x;
        input >> y;
        input >> z;

        if(input.eof())
            break;

        Vertex3D v = Vertex3D(x,y,z);
        mesh.addVertex(v);
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topsimplexes);

    //insert top simplexes
    for(unsigned long int i=0;i<num_topsimplexes;i++){
        int v[4];
        for(int j=0;j<4;j++)
            input >> v[j];
        if(v[0] != 3){
            cout<<"Where's the 3??? Not a valid .off file: "<<path<<endl;
            return false;
        }
        Triangle t = Triangle(v[1],v[2],v[3]);
        mesh.addTopSimplex(t);
    }
    return true;
}
