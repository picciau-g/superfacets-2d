#ifndef _READER_H
#define	_READER_H

#include <vector>
#include <string>

//#include "Vertex2D.h"
#include "Vertex3D.h"
#include "Mesh.h"
#include "Triangle.h"

using namespace std;

/// A class that read file and initializite some relevant library structures
class Reader {
public:
    /// A public method that reads a file containing a triangular mesh
    /*!
     * \param mesh a Mesh<T,V>& argument, represent the mesh to initialize
     * \param path a string argument, represent the path to the mesh file
     * \return a boolean value, true if the file is correctly readed, false otherwise
     */
    //static bool readMeshFile(Mesh<Vertex2D,Triangle>& mesh, string path);
    static bool readMeshFile(Mesh<Vertex3D,Triangle>& mesh, string path);
    static bool readOFFMesh(Mesh<Vertex3D, Triangle>& mesh, string path);
private:
    ///A constructor method
    Reader();
    ///A constructor method
    Reader(const Reader& orig);
    ///A destructor method
    virtual ~Reader();
};

#endif	/* _READER_H */

