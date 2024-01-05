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

#include "segmenter.h"
#include <stdio.h>


Segmenter::Segmenter(const Mesh<Triangle>& pMesh, unsigned int pRegions, float pAlpha)
    :
    m_Mesh(pMesh)
    , m_NCluster(pRegions)
    , m_Alpha(pAlpha)
    , m_FloodInit(false)
{
    std::cout << "Segmenter Given Number of Regions " << std::endl;
}

Segmenter::Segmenter(const Mesh<Triangle>& pMesh, float pRegionSize, float pAlpha, bool pFlood)
    :
    m_Mesh(pMesh)
    , m_Alpha(pAlpha)
    , m_RegionRadius(pRegionSize)
    , m_FloodInit(pFlood)
{
    std::cout << "Segmenter Given Expected Radius of Regions " << std::endl;
}


/**
 * @brief Segmenter::StartSegmentation loads the m_Mesh and initialize the structures
 */
void Segmenter::StartSegmentation()
{
    //initialize the structures
    m_FacesCentroids=ComputeFacesCentroids();
    m_NearestT=NearestFace();
    m_Normals.reserve(m_Mesh.GetNumberOfTopSimplexes()*sizeof(Normals));

    m_ClusterIndex = std::vector<int>(m_Mesh.GetNumberOfTopSimplexes(), -1);

    //Normals
    for(unsigned int a=0; a<m_Mesh.GetNumberOfTopSimplexes(); a++)
    {
        Triangle T=m_Mesh.GetTopSimplex(a);
        Normals n=Normals(m_Mesh.GetVertex(T.TV(0)),m_Mesh.GetVertex(T.TV(1)),m_Mesh.GetVertex(T.TV(2)));
        m_Normals.push_back(n);
    }

    //bounding box diagonal
    GetBoundingBoxDiagonal();

    BuildGlobalDistances();

    if(m_DebugMode)
    {
        std::cout << "start timer" << std::endl;
        std::cout<< "Clusters "<< m_NCluster << endl;
    }
    m_Timer.start();
    if(m_NCluster <= 0)
    { /// we don't know how many regions because it is radius-based

        if(m_FloodInit)
            FloodInitialization(m_NearestT);
        else
            GridInitialization();
    }
    else
    { /// We pass the number of regions as a parameter

        if(m_DebugMode)
            std::cout << "start seed placement" << std::endl;
        double meshArea = m_Mesh.MeshArea();
        double auxRad = sqrt(meshArea/(m_NCluster*M_PI));

        m_RegionRadius = auxRad/m_AABBDiagonal;
        PlaceSeeds(m_NearestT);
        if(m_DebugMode)
            std::cout << "MaxD " << m_RegionRadius << endl;
    }
    m_Timer.stop();

    initTime = m_Timer.getElapsedTimeInMilliSec();
    if(m_DebugMode)
        std::cout << "Time for the initialization is " << initTime << " milliseconds" <<std::endl;

    if(m_NCluster <= 0)
    { /// was not passed as parameter

        m_NCluster = 0;
        for(int ii=0; ii<m_Mesh.GetNumberOfTopSimplexes(); ii++)
        {

            if(m_NCluster <= m_ClusterIndex[ii])
                m_NCluster++;
        }
    }
    m_HasMoved=true;
    actual_iteration=0;

    std::cout << "There are " <<m_NCluster << " regions" << std::endl;
}

/**
 * @brief Segmenter::NearestFace
 * @return the index of the fac closest to mesh barycenter
 */
int Segmenter::NearestFace()
{

    float dist = glm::distance2(m_Mesh.GetBarycenter(), m_FacesCentroids.at(0).GetCoordinates());
    int nearestIndex = 0;

    for(unsigned int ii=0; ii<m_FacesCentroids.size(); ii++)
    {
        auto currentDistance = glm::distance2(m_Mesh.GetBarycenter(), m_FacesCentroids.at(ii).GetCoordinates());
        if(dist > currentDistance)
        {
            dist = currentDistance;
            nearestIndex=ii;
        }
    }
    return nearestIndex;
}

//A single step of the segmentation (interactive mode)
void Segmenter::Segmentation()
{

    //Any difference?
    m_HasMoved = UpdateCenters();

    if(m_HasMoved)
        ClassificationStep();
    else
        std::cout<< "Already converged" << std::endl;
}

/**
 * @brief Segmenter::ComputeFacesCentroids
 * @return a vector storing the centroid of each face
 */
vector<Vertex3D> Segmenter::ComputeFacesCentroids()
{
    std::vector<Vertex3D> res;

    for(int i=0; i<m_Mesh.GetNumberOfTopSimplexes(); i++)
    {
        Triangle T=m_Mesh.GetTopSimplex(i);

        glm::vec3 centroidCoords(0.0f, 0.0f, 0.0f);

        for(int a=0; a<T.getVerticesNum(); a++)
        {
            centroidCoords += m_Mesh.GetVertex(T.TV(a)).GetCoordinates();
        }

        centroidCoords /= 3.0f;

        res.push_back(Vertex3D(centroidCoords));
    }

    return res;
}


/**
 * @brief Segmenter::CalculateCentroidDistance
 * @param f1 index of the first triangle
 * @param f2 index of the second triangle
 * @return approximation of the geodesic distance between
 *      the two triangle barycenters
 */
float Segmenter::CalculateCentroidDistance(int f1, int f2)
{

    //first triangle
    Triangle t1 = m_Mesh.GetTopSimplex(f1);

    int indexT = -1;
    for(int ii=0;ii<3;ii++)
    {
        if(t1.TT(ii)==f2)
        {
            indexT=f2;
            break;
        }
    }

    //There's something wrong, f1 and f2 are not adjacent
    if(indexT < 0)
    {
        std::cout << "Not a valid face adjacency " << std::endl;
        exit(-1);
    }

    // Common edge and vertices
    Edge* shared = t1.TE((indexT)%3);
    auto indV0 = shared->EV(0);
    auto indV1 = shared->EV(1);
    Vertex3D halfP = m_Mesh.HalfPoint(indV0, indV1);

    //Corresponding centroids
    Vertex3D C1 = m_FacesCentroids.at(f1);
    Vertex3D C2 = m_FacesCentroids.at(f2);

    return C1.distance(halfP) + halfP.distance(C2);
}

/**
 * @brief Segmenter::BuildFaceDistances
 * @return a hash map with the approximate geodesic distances between
 *      pairs of adjacent triangles
 */
std::unordered_map<edgekey, float> Segmenter::BuildFaceDistances()
{

    std::unordered_map<edgekey, float> FD;

    for(unsigned int ii=0; ii<m_Mesh.GetNumberOfTopSimplexes(); ii++)
    {

        Triangle T = m_Mesh.GetTopSimplex(ii);

        for(int jj=0; jj<3; jj++)
        {

            int f2 = T.TT(jj);
            edgekey ek = KeyHandler::getKey(ii, f2);

            //If it is not already in the structure
            if(!FD.count(ek) && f2 >= 0)
            {
                float dist = CalculateCentroidDistance(ii, f2);
                assert(dist > 0);
                FD[ek] = dist;
            }
        }
    }
    return FD;
}

/**
 * @brief Segmenter::BuildAngleDistances
 * @return a hash map with the angular distance between
 *      pairs of adjacent triangles
 */
std::unordered_map<edgekey, float> Segmenter::BuildAngleDistances()
{

    std::unordered_map<edgekey, float> aDistances;

    for(int ii=0;ii<m_Mesh.GetNumberOfTopSimplexes();ii++)
    {
        Triangle T=m_Mesh.GetTopSimplex(ii);

        for(int jj=0;jj<3;jj++)
        {
            int f2 = T.TT(jj);
            edgekey ek = KeyHandler::getKey(ii,f2);

            if(aDistances.count(ek)==0 && f2>=0)
            {
                //Vertices of the shared edge
                faceind vi1 = T.TV((jj+1)%3);
                faceind vi2 = T.TV((jj+2)%3);

                //Length of shared edge between the two triangles
                Vertex3D v=m_Mesh.GetVertex(vi1);
                Vertex3D w=m_Mesh.GetVertex(vi2);
                Vertex3D vDifference(v.GetCoordinates()-w.GetCoordinates());
                auto balanceF = glm::length(vDifference.GetCoordinates());

                Normals N = m_Normals.at(ii);

                Vertex3D C1=m_FacesCentroids.at(ii);
                Vertex3D C2=m_FacesCentroids.at(f2);

                glm::vec3 diffCenters = C2.GetCoordinates() - C1.GetCoordinates();


                Normals nn;
                nn.SetNormal(diffCenters);


                float eta = (N.DotProd(nn) > 0) ? 1.0 : ETA_CONVEX;

                float dot = N.DotProd(m_Normals.at(f2));
                //Clamp dot product to be in [-1, 1]
                if(dot > 1.0)
                    dot=1.0;
                if(dot < -1.0)
                    dot=-1.0;

                //Divide it by pi to have it normalized in [0,1]
                float arcC=acos(dot)/M_PI;
                float d_a=eta*(arcC)*balanceF;

                aDistances[ek]=d_a;
            }
        }
    }
    return aDistances;
}


/**
 * @brief Segmenter::BuildGlobalDistances builds the hash map which stores the combined
 *          distance measure between pairs of adjacent triangles
 */
void Segmenter::BuildGlobalDistances()
{
    auto geodesicDistances = BuildFaceDistances();
    auto angularDistances = BuildAngleDistances();


    for(unsigned int ii=0; ii<m_Mesh.GetNumberOfTopSimplexes(); ii++)
    {

        Triangle T = m_Mesh.GetTopSimplex(ii);

        for(int jj=0; jj<3; jj++)
        {
            int f2=T.TT(jj);
            edgekey ek = KeyHandler::getKey(ii, f2);

            if(f2 >= 0 && !m_GlobalDistances.count(ek))
            {
                m_GlobalDistances[ek] = (geodesicDistances[ek] + m_Alpha*angularDistances[ek]) / m_AABBDiagonal;
            }
        }
    }

    geodesicDistances.clear();
    angularDistances.clear();
}

/**
 * @brief Segmenter::ClassificationStep performs a Dijkstra-based expansion
 */
void Segmenter::ClassificationStep()
{

    // initialize graph distance to infinity and cluster index to unassigned
    for(unsigned int ii=0;ii<m_FacesCentroids.size();ii++)
    {
        m_OutputDijkstra[ii]=FLT_MAX;
        m_ClusterIndex[ii] = -1;
    }

    for(auto regC : m_RegionCentroids) //initializa the region centroids
    {
        m_OutputDijkstra[regC.second]=0.0;    //  distance to centeroid is always 0
    }


    for(auto regC : m_RegionCentroids)//for(unsigned int ii=0;ii<m_RegionCentroids.size();ii++)
    {
        //Visit the triangles with a Dijkstra-based algorithm starting from a region centroid
        ExpandFromSeed(regC.second, regC.first);
    }
    std::cout << "Expanded " << std::endl;

    //Repair step to re-assign faces whose index is -1
    while(!CheckClusterIndex())
    {
        int violator = -1;
        for(size_t vv=0; vv<m_FacesCentroids.size(); vv++)
        {
            if(m_ClusterIndex[vv]<0)
            {
                violator=vv;
                break;
            }
        }

        if(violator >= 0 && violator < m_FacesCentroids.size())
        {
            ExpandFromSeed(violator, m_NCluster++);
        }
    }
    assert(CheckClusterIndex());
}

/**
 * @brief Segmenter::ExpandFromSeed expands a single centroid
 * @param indexT index of the triangle which is the center of the current region
 * @param newind index of the current region
 */
void Segmenter::ExpandFromSeed(int indexT, int newind)
{

    std::priority_queue<pointDist, std::vector<pointDist>, compare> Q;
    std::unordered_set<faceind> visited;

    pointDist seed;
    seed.indexP=indexT;
    seed.distanceP=0.0;
    Q.push(seed);

    m_ClusterIndex[indexT]=newind;

    m_RegionCentroids[newind]=indexT;

    while(!Q.empty())
    {

        pointDist current=Q.top();
        Q.pop();

        auto compNum = m_FacesCentroids[indexT].distance(m_FacesCentroids[current.indexP]);
        double comp = compNum/m_AABBDiagonal;
        if(comp <= 2.0*m_RegionRadius)
        {
            int neigh;
            Triangle T = m_Mesh.GetTopSimplex(current.indexP);

            for(int jj=0; jj<3; jj++)
            {
                neigh=T.TT(jj);

                if(neigh>= 0 && visited.count(neigh)==0)
                {
                    float newdist=current.distanceP + m_GlobalDistances[KeyHandler::getKey(current.indexP, neigh)];

                    if(newdist < m_OutputDijkstra[neigh])
                    {
                        m_OutputDijkstra[neigh]=newdist;
                        m_ClusterIndex[neigh]=newind;

                        pointDist pd;
                        pd.indexP=neigh;
                        pd.distanceP=newdist;
                        Q.push(pd);
                    }
                }
            }
            visited.insert(current.indexP);
        }
    }
}

/**
 * @brief Segmenter::PlaceSeeds Initialization if we pass the number of desired region as parameter
 * @param index index of the first centroid
 */
void Segmenter::PlaceSeeds(int index)
{

    m_RegionCentroids[0] = index;
    m_ClusterIndex[index]=0;

    //set with the face index of the already placed centroids
    std::unordered_set<faceind> aux_centroids;
    aux_centroids.insert(index);

    int count=1;
    if(m_DebugMode)
        std::cout << "Initializing..." << std::endl;

    while(m_RegionCentroids.size() < m_NCluster)
    {
        std::vector<double> distanceFromSeeds(m_Mesh.GetNumberOfTopSimplexes(), DBL_MAX);
        for(int iterFaces = 0; iterFaces < m_Mesh.GetNumberOfTopSimplexes(); iterFaces++)
        {
            if(aux_centroids.count(iterFaces) == 0)
            {
                distanceFromSeeds[iterFaces] = m_FacesCentroids[iterFaces].distance(m_FacesCentroids[m_RegionCentroids[0]]);

                for(size_t kk=1; kk<m_RegionCentroids.size(); kk++)
                {
                    double dd = m_FacesCentroids[iterFaces].distance(m_FacesCentroids[m_RegionCentroids[kk]]);
                    if(distanceFromSeeds[iterFaces] > dd)
                    {
                        distanceFromSeeds[iterFaces] = dd;
                    }

                }
            }
            else
            {
                distanceFromSeeds[iterFaces] = 0.0;
            }
        }

        faceind indOfM = IndexOfMax(distanceFromSeeds);
        m_ClusterIndex[indOfM] = count;
        aux_centroids.insert(indOfM);
        m_RegionCentroids[count++] = indOfM;
    }

    if(m_DebugMode)
        std::cout << "Done: regions " << m_RegionCentroids.size() << std::endl;
}

/**
 * @brief Segmenter::GridInitialization initialize the segmentation dividing the object with a regular grid
 */
void Segmenter::GridInitialization()
{

    std::unordered_map<gridkey, int> hmap; //for already assigned maps
    int id=0; //index of current region
    float denominator = 2*m_RegionRadius*m_AABBDiagonal;

    for(unsigned int ii=0;ii<m_FacesCentroids.size();ii++)
    {

        Vertex3D currentV = m_FacesCentroids[ii];
        //find which region it belongs to
        glm::vec3 gridCell = currentV.GetCoordinates() / denominator;

        //if it is the first time we encounter the region
        if(hmap.count(KeyHandler::getGridKey(gridCell.x,gridCell.y,gridCell.z))==0)
        {
            m_RegionCentroids[id]=ii;
            hmap[KeyHandler::getGridKey(gridCell.x,gridCell.y,gridCell.z)] = id++;
        }

        //assign a region index to the triangle
        m_ClusterIndex[ii]=hmap[KeyHandler::getGridKey(gridCell.x,gridCell.y,gridCell.z)];
    }
}

/**
 * @brief Segmenter::FloodInitialization initialization with a Dijkstra-based expansion
 * @param indexT index of the startin point
 */
void Segmenter::FloodInitialization(int indexT)
{
    priority_queue<pointDist, vector<pointDist>, compare> facesQueue;

    // Initialize the distance graph to infinity for each triangle
    for(unsigned int ii=0; ii<m_Mesh.GetNumberOfTopSimplexes(); ii++)
    {
        m_OutputDijkstra[ii]=FLT_MAX;
    }
    m_OutputDijkstra[indexT]=0; //Distance of the seed is 0

    //Queue that stores all the triangles ordered with respect to their Euclidean distance from the seed
    //we had it to avoid having unassigned faces at the end of the initialization
    for(unsigned int a=0; a<m_Mesh.GetNumberOfTopSimplexes(); a++)
    {
        pointDist pDist;
        pDist.indexP=a;
        if(a == indexT)
            pDist.distanceP=0.0;
        else
            pDist.distanceP=m_FacesCentroids[a].distance(m_FacesCentroids[indexT]);

        facesQueue.push(pDist);
    }


    //Index of the current region
    int currentIndex=0;
    // Main loop
    while(!facesQueue.empty())
    {

        // Remove faces already assigned to some cluster
        while(facesQueue.size() >= 1 && m_ClusterIndex[facesQueue.top().indexP]>=0)
            facesQueue.pop();

        //Take the first element of the queue
        pointDist center=facesQueue.top();
        facesQueue.pop();

        if(center.distanceP==FLT_MAX)
        { //Unconnected or unvisited face
            break;
        }

        ExpandFromSeed(center.indexP, currentIndex++);
    }
}

/**
 * @brief Segmenter::GetBoundingBoxDiagonal sets lower left and upper right corners of the m_Mesh bounding box
 *          and its diagonal length
 */
void Segmenter::GetBoundingBoxDiagonal()
{
    m_AABBDiagonal = m_Mesh.GetGoundingBoxDiagonal();
}

/**
 * @brief Segmenter::CheckClusterIndex
 * @return true if every triangle is assigned to some cluster (its index is >=0)
 */
bool Segmenter::CheckClusterIndex()
{
    for(int a=0; a<m_Mesh.GetNumberOfTopSimplexes(); a++)
    {
        if(m_ClusterIndex[a] < 0)
        {
            return false;
        }
    }

    return true;
}


/**
 * @brief Segmenter::UpdateCenters takes centroids as close as possible to the center of the region
 * @return true if no centroid has moved from the previous iteration
 */
bool Segmenter::UpdateCenters()
{

    // differences in each region
    std::vector<int >differences(m_NCluster, 0);

    bool any_moves=false;

    if(m_DebugMode)
    {
        std::cout << "Update" << endl;
    }

    // Iterate on each one of the regions
    for(unsigned int ii=0; ii<m_RegionCentroids.size(); ii++)
    {
        //Area of a region
        double regionArea=0.0;
        for(unsigned int ff=0;ff<m_FacesCentroids.size();ff++)
        {
            if(m_ClusterIndex[ff]==ii)
                regionArea += m_Mesh.TriangleArea(ff);
        }

        glm::vec3 currentCentroid = glm::vec3(0.0, 0.0, 0.0);

        // Area-weighted average
        for(unsigned int ff=0;ff<m_FacesCentroids.size();ff++)
        {
            if(m_ClusterIndex[ff] == ii)
            {
                glm::vec3 currentFace=m_FacesCentroids.at(ff).GetCoordinates() * static_cast<float>(m_Mesh.TriangleArea(ff));
                currentCentroid += currentFace;

            }
        }

        //Average of the coordinates of the baricenters of the region
        Vertex3D nc(currentCentroid / static_cast<float>(regionArea));

        double actualD=nc.distance(m_FacesCentroids[m_RegionCentroids[ii]]);

        //Centroid closest to nc
        for(unsigned int ff=0; ff<m_FacesCentroids.size(); ff++)
        {
            if(m_ClusterIndex[ff]==ii)
            {
                if(nc.distance(m_FacesCentroids.at(ff)) < actualD && ff != m_RegionCentroids[ii])
                {
                    actualD=nc.distance(m_FacesCentroids[ff]);
                    any_moves=true;
                    m_RegionCentroids[ii]=ff;
                    differences[ii]++;
                }
            }
        }
    }

    // count the differences
    int ctrdiff=0;
    int totaldiff=0;

    for(int k=0;k<m_NCluster;k++)
    {
        if(differences[k])
            ctrdiff++;
        totaldiff+=differences[k];
    }
    if(m_DebugMode)
    {
        std::cout << "Regions " << m_NCluster << std::endl;
        std::cout << "Differences in " << ctrdiff << " regions, total = " << totaldiff << std::endl;
    }
    return any_moves;
}


bool Segmenter::CheckIfAllVisited(std::unordered_set<edgekey> pSet)
{

    for(auto it=pSet.begin(); it!=pSet.end(); ++it)
    {
        edgekey faceV = *it;
        if(faceV > m_ClusterIndex.size() || m_ClusterIndex[faceV] < 0)
            return false;
    }
    return true;
}

faceind Segmenter::IndexOfMax(const std::vector<double> &pArray)
{
    double current = DBL_MIN; //All the values should be > 0 (they are distances)
    int toRet=-1;
    for(int aa=0; aa<pArray.size(); aa++)
    {
        if(pArray[aa] > current)
        {
            toRet = aa;
            current = pArray[aa];
        }
    }
    return toRet;
}
