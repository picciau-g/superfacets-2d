#ifndef VERTEXBASEDSEGMENTER_H
#define VERTEXBASEDSEGMENTER_H

//#include "normals.h"
//#include "Timer.h"
//#include "Reader.h"
//#include <tr1/unordered_map>
//#include <queue>
//#include <tr1/unordered_set>
//#include <numeric>
//#include <QString>
//#include <QStringList>
//#include <cfloat>
//#include <fstream>
#include "common.h"

//namespace std { using namespace __gnu_cxx; }

typedef unsigned long int vertexind;
//typedef unsigned long long int edgekey;

//struct pointDist{
//    int indexP;
//    double distanceP;
//};

//struct compare{

//    bool operator ()(pointDist p1, pointDist p2){

//        return p1.distanceP > p2.distanceP;
//    }
//};

edgekey getVertexKey(vertexind a, vertexind b);

class VertexBasedSegmenter
{
public:
    VertexBasedSegmenter();
    string filename;
    string fieldfilename;

    Mesh<Vertex3D, Triangle> mesh;
    vector<Normals> norms;
    vector<float> functionValue;

    std::unordered_map<edgekey,float> vertexDistances;
    std::unordered_map<edgekey,float> functionVDistances;
    std::unordered_map<edgekey,float> globalDistances;

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

    inline void setAlpha(float a){
        this->alpha=a;
    }

    inline void startSeg(){
        this->Segmentation();
    }

    inline void callWriter(string outfile){
        writeSegmOnFile(outfile);
    }

private:
    float *faceAreas;
    int *clusterIndex;
    int NCluster;
    float maxD;
    float alpha;

    /// Center of the mesh (it may not be on the mesh surface)
    Vertex3D centerMesh;
    /// Lower left and upper right vertices of the mesh bounding box
    Vertex3D minCoords;
    Vertex3D maxCoords;
    /// Priority queue used for Dijkstra (stores the distances from the center)
    priority_queue<pointDist, vector<pointDist>, compare> globalQ;

    /// Length of the bouding box diagonal
    float BBDiagonal;

    /// Index of the "Central Vertex" of the region
    std::unordered_map<edgekey, int> regionCentroids;

    std::unordered_map<edgekey, float> outputDijkstra;



    void loadMesh();
    Vertex3D centerCoordinate();
    int nearestVertex(); //closest to mesh barycenter

    /// Functions to build geodesic (face), angular and global distances
    std::unordered_map<edgekey, float> buildVertexDistances();
    std::unordered_map<edgekey, float> buildFunctionVDistances();
    void buildGlobalDistances();

    void placeSeeds(int index);

    void expansionStep();
    void expandSeed(int indexV, int newInd);

    /// Input/Output Functions
    int writeSegmOnFile(string);
    void openMeshFile(string);
    void openCurvatureFile(string);

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

    inline bool CheckClusterIndex(){
        bool ret=true;

        for(int a=0;a<mesh.getNumVertex();a++){
            if(clusterIndex[a] < 0){
                ret=false;
            }
        }
        return ret;
    }

    inline vertexind indexOfMax(double *array){
        double actual = -10000; //All the values should be > 0 (they are distances)
        int toRet=-1;
        for(int aa=0;aa<mesh.getNumVertex();aa++){
            if(array[aa] > actual){
                toRet = aa;
                actual = array[aa];
            }
        }
        return toRet;
    }

};

#endif // VERTEXBASEDSEGMENTER_H
