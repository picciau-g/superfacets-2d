/*
 *
 *   2014
 *   Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
 *   Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
 *                 Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
 *
 *   Title:          Fast and scalable mesh superfacets
 *   Submission to Pacific Graphics 2014
 *
 *
 **/
#ifndef SEGMENTER_H
#define SEGMENTER_H

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

namespace std { using namespace __gnu_cxx; }

///// types for indices of triangles and for adjacencies between them
typedef unsigned long int faceind;
//typedef unsigned long long int edgekey;

///// types for grid initialization
//typedef signed char coordind;
//typedef int gridkey;

///*!
// * \brief The pointDist struct: used for the priority queue in the Dijkstra-based expansion
// */
//struct pointDist{
//    int indexP;
//    double distanceP;
//};

///*!
// * \brief The compare struct: redefinition of the > operator to be able to use a priority queue
// */
//struct compare{

//    bool operator ()(pointDist p1, pointDist p2){

//        return p1.distanceP > p2.distanceP;
//    }
//};

/**
 *
 * @brief The Segmenter class decompose an object into small patches according to
 * a distance metric that is a qeighted sum of spatial and angular information
 */
class Segmenter
{
public:
    Segmenter();
    string filename;
    Mesh<Vertex3D, Triangle> mesh;

    vector<Vertex3D> facesCentroids;
    vector<Normals> norms;

    /// if we want it to be verbose (debug purposes)
    bool debugMode;

    /// ONLY FOR THE COUNT OF EXPANSIONS
    int actual_iteration;

    /// Total average running time
    double runningTime;

    /// Time to initialize the segmentation
    double initTime;

    /// Whether we are using flood initialization or the squares one
    bool floodInit;

    /// To see if it needs to print the header on the  output file
    bool putHeader;
    /// Map containing the geodesic distances between pairs of faces
    std::tr1::unordered_map<edgekey,float> faceDistances;
    /// Map containing the angular distances between pairs of faces
    std::tr1::unordered_map<edgekey, float> angleDistances;
    /// Map containing the global distances between pairs of faces (obtained combining geodesic and angular distance)
    std::tr1::unordered_map<edgekey, float> globalDistances;

    std::tr1::unordered_map<edgekey, float> outputDijkstra;

    /// Index of the "Central Triangle" of the region
    std::tr1::unordered_map<edgekey, int> regionCentroids;
    /// Areas ot the faces of the mesh
    float* faceAreas;

    /// The index of the cluster to which a face belongs
    int* clusterIndex;

    /// How many times each face is epxanded
    int* expanded;

    /// Sets the actual dimension of the search region in the expansion step
    inline void setTimesR(int t){
        this->timesR=t;
    }

    /// Retrieve the number of clusters (if we did not set it as a parameter)
    inline int getNCluster(){
        return this->NCluster;
    }

    /// Set number of desired regions (if we already know it)
    inline void setNCluster(int num){
        this->NCluster = num;
    }

    /*!
         * \brief setAlpha
         * \param aa the desired value of alpha
         */
    inline void setAlpha(float aa){
        this->alpha=aa;
    }

    /*!
         * \brief setMaxD
         * \param dd the desired value for the search region radius
         */
    inline void setMaxD(float dd){
        this->maxD=dd;
    }

    inline float getMaxD(){
        return this->maxD;
    }

    /*!
         * \brief setEtaConvex
         * \param eta the desired value of eta (only if the angle is convex)
         */
    inline void setEtaConvex(double eta){
        this->etaconvex=eta;
    }

    /**
     * @brief callLoad wrapper to load the mesh from outside the class
     */
    inline void callLoad(){
        this->loadMesh();
    }

    /**
     * @brief setMaxIters
     * @param i Maximum number of desired iterations (default=50)
     */
    inline void setMaxIters(int i){
        this->maxIters = i;
    }

    inline int getMaxIters(){
        return this->maxIters;
    }

    /*!
         * \brief centerCoordinate
         * \return The center of the mesh (it can be outside the mesh surface)
         */
    Vertex3D centerCoordinate();

    /*!
         * \return the index of the nearest face to the center of the mesh
         */
    int nearestFace();

    /// Functions to build geodesic (face), angular and global distances
    std::tr1::unordered_map<edgekey, float> buildFaceDistances();
    std::tr1::unordered_map<edgekey, float> buildAngleDistances();
    void buildGlobalDistances();

    /// Alternative initialization that takes as input the number of regions
    void placeSeeds(int index);

    /// Dijkstra functions (initialization and iterative step)
    void initializationGrid();
    void floodInitialization(int);
    void expansionStep();
    void expandSeed(int, int);

    /// Input/Output Functions
    int writeSegmOnFile(string);

    /**
     * @brief setElapsedTime (for the .seg header)
     * @param t execution time in milliseconds
     */
    inline void setElapsedTime(double t){
        this->millisecs=t;
    }

    /**
     * @brief setIters (for the .seg header)
     * @param count number of iterations
     */
    inline void setIters(int count){
        this->iters=count;
    }

    /**
     * @brief setTypeVis
     * @param v
     *TODO
     */
    inline void setTypeVis(unsigned int v){
        this->visT=v;
    }

    /// Check that all the faces have been assigned to some cluster
    bool CheckClusterIndex();
    /// Update of the centroids
    bool updateCenters();

private:
    Timer TM;

    /// Time needed to perform the segmentation
    double millisecs;
    /// Number of iteration before the algorithm converges
    int iters;

    /// Center of the mesh (it may not be on the mesh surface)
    Vertex3D centerMesh;
    /// Lower left and upper right vertices of the mesh bounding box
    Vertex3D minCoords;
    Vertex3D maxCoords;

    /// Index of the face closest to the center
    int nearestT;
    /// Maximum distance within the same region
    float maxDistance;
    /// Number of regions
    int NCluster;

    /// How much angular distance will be important
    double alpha;
    /// The radius of the search region
    float maxD;

    /// Priority queue used for Dijkstra (stores the distances from the center)
    priority_queue<pointDist, vector<pointDist>, compare> globalQ;

    /// Length of the bouding box diagonal
    float BBDiagonal;

    /// Loads the mesh and initializes the structures
    void loadMesh();

    /// Performs a step of the segmentation
    void Segmentation();

    /// Gets the length of the bounding box
    void getBBDiagonal();

    /// Initializes the queue for Dijkstra algorithm (all distance to infinite except the starting one)
    priority_queue<pointDist, vector<pointDist>, compare> initializeQueue(int);

    inline int retFirst(int *array){
        //int toret = -1;
        for(int ii=0;ii<mesh.getTopSimplexesNum();ii++){
            if(array[ii]==0)
                return ii;
        }
        return -1;
    }

    /*!
     * \brief getKey
     * \param a index of the first face
     * \param b index of the second face
     * \return the key value (used by the structures storing the distances between faces)
     */
    inline edgekey getKey(faceind a, faceind b){

        if(a<=b)
            return (edgekey(a) << 32 | edgekey(b));
        else
            return (edgekey(b) << 32 | edgekey(a));
    }

    inline gridkey getGridKey(coordind x, coordind y, coordind z){
        //return (gridkey(x) << 16 | gridkey(y) << 8 | gridkey(z));
        if(x <= y && x <= z){
            if(y<=z)
                return (gridkey(fabs(x)) << 16 | gridkey(fabs(y)) << 8 | gridkey(fabs(z)));
            else
                return (gridkey(fabs(x)) << 16 | gridkey(fabs(z)) << 8 | gridkey(fabs(y)));
        }
        else{
            if(y <= z){
                if(x <= z)
                    return (gridkey(fabs(y)) << 16 | gridkey(fabs(x)) << 8 | gridkey(fabs(z)));
                return (gridkey(fabs(y)) << 16 | gridkey(fabs(z)) << 8 | gridkey(fabs(x)));
            }
            if(x <= y)
                return (gridkey(fabs(z)) << 16 | gridkey(fabs(x)) << 8 | gridkey(fabs(y)));
            return (gridkey(fabs(z)) << 16 | gridkey(fabs(y)) << 8 | gridkey(fabs(x)));
        }

    }

    /*!
         * \brief setAreas
         * Builds the array containing the area of the faces
         */
    inline void setAreas(){
        for(int a=0;a<mesh.getTopSimplexesNum();a++)
            faceAreas[a]=mesh.TArea(a);
    }

    /**
     * @brief indexOfMax
     * @param array array of double values
     * @return the index of the highest value in array
     */
    inline faceind indexOfMax(double *array){
        double actual = -10000; //All the values should be > 0 (they are distances)
        int toRet=-1;
        for(int aa=0;aa<mesh.getTopSimplexesNum();aa++){
            if(array[aa] > actual){
                toRet = aa;
                actual = array[aa];
            }
        }
        return toRet;
    }

    /**
     * @brief checkVisited check if all the faces were assigned
     * @param set set of faces
     * @return true if all faces were assigned, false otherwise
     */
    inline bool checkVisited(std::tr1::unordered_set<edgekey> set){

        for(auto it=set.begin(); it!=set.end(); ++it){
            edgekey faceV = *it;
            if(clusterIndex[faceV] < 0)
                return false;
        }
        return true;
    }

    /*!
         * \brief computeCentroids
         * \return A vector with all the centroids of the mesh faces
         */
    vector<Vertex3D> computeCentroids();

    /*!
     * \brief openMeshFile reads the extension of the mesh file (.off or .tri) and loads it
     */
    void openMeshFile(string);

    /*!
         * \brief halfPoint
         * \return The middle point of a segment
         */
    Vertex3D halfPoint(Vertex3D, Vertex3D);

    /// Returns the distance between two triangles (between their centers)
    float centroidDistance(int, int);

    /// value of eta if the angle is convex
    double etaconvex;

    /// true if any of the centroids has moved
    bool moves;

    /// The kind of visualization (normal, highliting concavities or fading)
    unsigned int visT;

    /// This value multiplies the search radius in the expansion step
    int timesR;

    /// iterative steps to be performed at most
    int maxIters;
};

#endif // SEGMENTER_H
