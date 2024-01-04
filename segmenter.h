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

#include "common.h"

#define ETA_CONVEX 2

class KeyHandler
{
public:
/*!
     * \brief getKey
     * \param a index of the first face
     * \param b index of the second face
     * \return the key value (used by the structures storing the distances between faces)
     */
static edgekey getKey(faceind a, faceind b)
{

    if(a<=b)
        return (edgekey(a) << 32 | edgekey(b));
    else
        return (edgekey(b) << 32 | edgekey(a));
}


static gridkey getGridKey(coordind x, coordind y, coordind z)
{
    //return (gridkey(x) << 16 | gridkey(y) << 8 | gridkey(z));
    if(x <= y && x <= z)
    {
        if(y<=z)
            return (gridkey(fabs(x)) << 16 | gridkey(fabs(y)) << 8 | gridkey(fabs(z)));
        else
            return (gridkey(fabs(x)) << 16 | gridkey(fabs(z)) << 8 | gridkey(fabs(y)));
    }
    else
    {
        if(y <= z)
        {
            if(x <= z)
                return (gridkey(fabs(y)) << 16 | gridkey(fabs(x)) << 8 | gridkey(fabs(z)));
            return (gridkey(fabs(y)) << 16 | gridkey(fabs(z)) << 8 | gridkey(fabs(x)));
        }
        if(x <= y)
            return (gridkey(fabs(z)) << 16 | gridkey(fabs(x)) << 8 | gridkey(fabs(y)));
        return (gridkey(fabs(z)) << 16 | gridkey(fabs(y)) << 8 | gridkey(fabs(x)));
    }

}

static int minArrayIndex(float array[3])
{
    if(array[0]<array[1] && array[0]<array[2])
        return 0;
    else if(array[1]<array[2])
        return 1;
    else
        return 2;
}

};



/**
 *
 * @brief The Segmenter class decompose an object into small patches according to
 * a distance metric that is a qeighted sum of spatial and angular information
 */
class Segmenter
{
public:
    Segmenter() = delete;
    Segmenter(const Segmenter& ) = delete;

    Segmenter(const Mesh<Triangle>& pMesh, unsigned int pRegions, float pAlpha, bool pFlood);
    Segmenter(const Mesh<Triangle>& pMesh, float pRegionSize, float pAlpha);

    virtual ~Segmenter(){}


    inline Mesh<Triangle> GetMesh() const
    {
        return m_Mesh;
    }

    inline std::string GetMeshName() const
    {
        return m_MeshName;
    }

    inline void SetMeshName(const std::string& pMeshName)
    {
        m_MeshName = pMeshName;
    }


    inline void SetPutHeader(bool pH)
    {
        m_PutHeader = pH;
    }

    /// Initializes the structures and starts the segmentation
    void StartSegmentation();
    void ClassificationStep();
    /// Update of the centroids
    bool UpdateCenters();

    /// Input/Output Functions
    int WriteSegmOnFile(string);


private:
    /*!
         * \return the index of the nearest face to the center of the mesh
         */
    int NearestFace();

    /// Functions to build geodesic (face), angular and global distances
    std::unordered_map<edgekey, float> BuildFaceDistances();
    std::unordered_map<edgekey, float> BuildAngleDistances();
    void BuildGlobalDistances();

    /// Alternative initialization that takes as input the number of regions
    void PlaceSeeds(int index);

    /// Dijkstra functions (initialization and iterative step)
    void GridInitialization();
    void FloodInitialization(int);
    void ExpandFromSeed(int, int);

    /// Performs a classification step
    void Segmentation();
    /// Check that all the faces have been assigned to some cluster
    bool CheckClusterIndex();

    /// Gets the length of the bounding box
    void GetBoundingBoxDiagonal();

    vector<Vertex3D> ComputeFacesCentroids();
    //MOVE IT TO WINDOW CLASS
    void OpenMeshFile(string);
    /// Returns the distance between two triangles (between their centers)
    float CalculateCentroidDistance(int, int);

    bool CheckIfAllVisited(std::unordered_set<edgekey> pSet);
    faceind IndexOfMax(const std::vector<double>& pArray);

    /// Initializes the queue for Dijkstra algorithm (all distance to infinite except the starting one)
    priority_queue<pointDist, vector<pointDist>, compare> InitializeDijkstraQueue(int);


private:
    std::string m_MeshName;
    Mesh<Triangle> m_Mesh;
    /// Map containing the global distances between pairs of faces (obtained combining geodesic and angular distance)
    std::unordered_map<edgekey, float> m_GlobalDistances;
    /// Distance between a face and its closest region centroid
    std::unordered_map<edgekey, float> m_OutputDijkstra;
    /// Region index for each face THIS WILL LIKELY BE DELETED AND MERGED IN OUTPUTDIJKSTRA
    std::vector<int> m_ClusterIndex;

    /// Index of the "Central Triangle" of the region
    std::map<edgekey, int> m_RegionCentroids;

    std::vector<Vertex3D> m_FacesCentroids;
    std::vector<Normals> m_Normals;
    Timer m_Timer;

    /// if we want it to be verbose (debug purposes)
    bool m_DebugMode;

    /// ONLY FOR THE COUNT OF EXPANSIONS
    int actual_iteration;

    /// Total average running time
    double runningTime;

    /// Time to initialize the segmentation
    double initTime;

    /// Index of the face closest to the center
    int m_NearestT;
    /// Number of regions
    int m_NCluster;

    /// How much angular distance will be important
    double m_Alpha;
    /// The radius of the search region
    float m_RegionRadius;

    /// Length of the bouding box diagonal
    float m_AABBDiagonal;

    /// true if any of the centroids has moved
    bool m_HasMoved;
    /// Whether we are using flood initialization or the squares one
    bool m_FloodInit;
    /// To see if it needs to print the header on the  output file
    bool m_PutHeader;
};

#endif // SEGMENTER_H
