#ifndef VERTEXBASEDSEGMENTER_H
#define VERTEXBASEDSEGMENTER_H

#include "normals.h"
#include "Timer.h"
#include "Reader.h"
#include <tr1/unordered_map>
#include <queue>
#include <tr1/unordered_set>
#include <numeric>
#include <QString>
#include <QStringList>
#include <cfloat>
#include <fstream>

namespace std { using namespace __gnu_cxx; }

typedef unsigned long int vertexind;
typedef unsigned long long int edgekey;

struct pointDist{
    int indexP;
    double distanceP;
};

struct compare{

    bool operator ()(pointDist p1, pointDist p2){

        return p1.distanceP > p2.distanceP;
    }
};

edgekey getKey(faceind a, faceind b){

        if(a<=b)
            return (edgekey(a) << 32 | edgekey(b));
        else
            return (edgekey(b) << 32 | edgekey(a));
    }

class VertexBasedSegmenter
{
public:
    VertexBasedSegmenter();
    string filename;

    Mesh<Vertex3D, Triangle> mesh;
    vector<Normals> norms;
    vector<float> functionValue;

    std::tr1::unordered_map<edgekey,float> vertexDistances;
    std::tr1::unordered_map<edgekey,float> functionVDistances;
    std::tr1::unordered_map<edgekey,float> globalDistances;

    inline int getNCluster(){
        return this->NCluster;
    }
    inline void setNCluster(int NC){
        this->NCluster = NC;
    }

    inline void setMaxD(float dd){
        this->maxD=dd;
    }

    inline float getMaxD(){
        return this->maxD;
    }

    inline void callLoad(){
        this->loadMesh();
    }

private:
    float *faceAreas;
    int *clusterIndex;
    int NCluster;
    float maxD;

    /// Center of the mesh (it may not be on the mesh surface)
    Vertex3D centerMesh;
    /// Lower left and upper right vertices of the mesh bounding box
    Vertex3D minCoords;
    Vertex3D maxCoords;
    /// Priority queue used for Dijkstra (stores the distances from the center)
    priority_queue<pointDist, vector<pointDist>, compare> globalQ;

    /// Length of the bouding box diagonal
    float BBDiagonal;


    void loadMesh();
    Vertex3D centerCoordinate();
    int nearestVertex(); //closest to mesh barycenter

    /// Functions to build geodesic (face), angular and global distances
    std::tr1::unordered_map<edgekey, float> buildVertexDistances();
    std::tr1::unordered_map<edgekey, float> buildFunctionVDistances();
    void buildGlobalDistances();

    void placeSeeds(int index);
    void startSeg();

    void expansionStep();
    void expandSeed(int, int);

    /// Input/Output Functions
    int writeSegmOnFile(string);
    void openMeshFile(string);

    /// Check that all the vertices have been assigned to some cluster
    bool CheckClusterIndex();
    /// Update of the centroids
    bool updateCenters();

    /// Performs a step of the segmentation
    void Segmentation();

    /// Gets the length of the bounding box
    void getBBDiagonal();

    /// Initializes the queue for Dijkstra algorithm (all distance to infinite except the starting one)
    priority_queue<pointDist, vector<pointDist>, compare> initializeQueue(int);

    inline void setAreas(){
        for(int a=0;a<mesh.getTopSimplexesNum();a++)
            faceAreas[a]=mesh.TArea(a);
    }
    Vertex3D halfPoint(Vertex3D, Vertex3D);

};

#endif // VERTEXBASEDSEGMENTER_H
