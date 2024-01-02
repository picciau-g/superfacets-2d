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

/**
 * @brief Segmenter::Segmenter Constructor
 */
Segmenter::Segmenter()
{
    //for running times
    TM=Timer();
    //initially it is 0
    runningTime=0.0;

    //If we don't put this, flood and grid do not work
    NCluster = -1;

    debugMode=false;
}

/**
 * @brief Segmenter::openMeshFile opens the file which stores the mesh and load the structure
 * @param mName path to the mesh file
 */
void Segmenter::openMeshFile(string mName){

    mesh = Mesh<Triangle>();
    QString qmn = QString::fromStdString(mName);
    QStringList qsl = qmn.split(".");

    if(!qsl.back().compare("tri"))
        Reader::readMeshFile(mesh, mName);
    else if(!qsl.back().compare("off"))
        Reader::readOFFMesh(mesh, mName);
    else{
        cout<<"Not a valid file format (it must be .off or .tri)"<<endl;
        exit(0);
    }

    mesh.Build();
}

/**
 * @brief Segmenter::loadMesh loads the mesh and initialize the structures
 */
void Segmenter::loadMesh(){

    TM.start();
    openMeshFile(filename);

    if(debugMode)
        cout<<"Time to load and build is "<<TM.getElapsedTimeInMilliSec()<<" milleseconds"<<endl;

    //initialize the structures
    facesCentroids=computeCentroids();
    centerMesh=centerCoordinate();
    nearestT=nearestFace();
    faceAreas=new float[mesh.GetNumberOfTopSimplexes()];
    norms.reserve(mesh.GetNumberOfTopSimplexes()*sizeof(Normals));

    //Normals
    for(unsigned int a=0;a<mesh.GetNumberOfTopSimplexes();a++){
        Triangle T=mesh.GetTopSimplex(a);
        Normals n=Normals(mesh.GetVertex(T.TV(0)),mesh.GetVertex(T.TV(1)),mesh.GetVertex(T.TV(2)));
        norms.push_back(n);
    }
    //Area of the faces
    setAreas();

    //bounding box diagonal
    getBBDiagonal();

    //geodesic, angular and global distances
    faceDistances = buildFaceDistances();
    angleDistances = buildAngleDistances();

    buildGlobalDistances();
    faceDistances.erase(faceDistances.begin(), faceDistances.end());
    angleDistances.erase(angleDistances.begin(), angleDistances.end());

    clusterIndex = new int[mesh.GetNumberOfTopSimplexes()];
    if(debugMode)
        expanded = new int[mesh.GetNumberOfTopSimplexes()]; // it will be used only for debug purposes

    //Initialization (all the three methods are possible)
    if(debugMode){
        cout<<"start timer"<<endl;
        cout<<"Clusters "<<NCluster<<endl;
    }
    TM.start();
    if(NCluster < 0){ /// we don't know how many regions because it is radius-based

        cout<<"Flood! "<<endl;
        if(floodInit)
            floodInitialization(nearestT);
        else
            initializationGrid();
    }
    else{ /// We pass the number of regions as a parameter

        cout<<"Seed"<<endl;
        if(debugMode)
            cout<<"start seed placement"<<endl;
        double meshArea = mesh.MeshArea();
        double auxRad = sqrt(meshArea/(NCluster*M_PI));

        this->maxD = auxRad/BBDiagonal;
        placeSeeds(nearestT);
        if(debugMode)
            cout<<"MaxD "<<maxD<<endl;
    }
    TM.stop();

    initTime = TM.getElapsedTimeInMilliSec();
    if(debugMode)
        cout<<"Time for the initialization is "<<initTime<<" milliseconds"<<endl;

    if(NCluster < 0){ /// was not passed as parameter

        NCluster = 0;
        for(int ii=0; ii<mesh.GetNumberOfTopSimplexes(); ii++){

            if(NCluster <= clusterIndex[ii])
                NCluster++;
        }
    }
    moves=true;
    actual_iteration=0;

    cout<<"There are "<<NCluster<<" regions"<<endl;
}

/**
 * @brief Segmenter::centerCoordinate
 * @return return the x,y,and z coordinates of the mesh barycenter
 */
Vertex3D Segmenter::centerCoordinate(){

    Vertex3D res;
    float sumX=0.0; float sumY=0.0; float sumZ=0.0;

    for(int ii=0;ii<mesh.GetNumVertices();ii++){

        sumX += mesh.GetVertex(ii).getX();
        sumY += mesh.GetVertex(ii).getY();
        sumZ += mesh.GetVertex(ii).getZ();
    }
    res.setX(sumX/mesh.GetNumVertices());
    res.setY(sumY/mesh.GetNumVertices());
    res.setZ(sumZ/mesh.GetNumVertices());
    if(debugMode)
        cout<<"RES "<<res.getX()<<" "<<res.getY()<<" "<<res.getZ()<<endl;
    return res;
}

/**
 * @brief Segmenter::nearestFace
 * @return the index of the fac closest to mesh barycenter
 */
int Segmenter::nearestFace(){

    float dist = centerMesh.distance(facesCentroids.at(0));
    int nearestIndex;

    for(unsigned int ii=0; ii<facesCentroids.size(); ii++){
        if(dist > centerMesh.distance(facesCentroids.at(ii))){
            dist = centerMesh.distance(facesCentroids.at(ii));
            nearestIndex=ii;
        }
    }
    return nearestIndex;
}

//A single step of the segmentation (interactive mode)
void Segmenter::Segmentation(){

    //Any difference?
    moves = updateCenters();

    if(moves)
        expansionStep();
    else
        cout<<"Already converged"<<endl;
}

/**
 * @brief Segmenter::computeCentroids
 * @return a vector storing the centroid of each face
 */
vector<Vertex3D> Segmenter::computeCentroids(){

    vector<Vertex3D> res;
    Vertex3D aux;

    for(int i=0;i<mesh.GetNumberOfTopSimplexes();i++){
        Triangle T=mesh.GetTopSimplex(i);
        int indices[3];
        for(int a=0;a<3;a++)
            indices[a]=T.TV(a);
        Vertex3D v3d[3];
        for(int a=0;a<3;a++){
            v3d[a]=mesh.GetVertex(indices[a]);
        }

        float valX=0.0, valY=0.0, valZ=0.0;
        for(int a=0;a<3;a++){
            valX+=v3d[a].getX();
            valY+=v3d[a].getY();
            valZ+=v3d[a].getZ();
        }
        aux.setX(valX/3);
        aux.setY(valY/3);
        aux.setZ(valZ/3);
        res.push_back(aux);
    }
    return res;
}

/**
 * @brief Segmenter::halfPoint half point of an edge
 * @param v1 first vertex
 * @param v2 second vertex
 * @return the half point of edge v1v2
 */
Vertex3D Segmenter::halfPoint(Vertex3D v1, Vertex3D v2){

    double coordX = (v1.getX() + v2.getX())/2;
    double coordY = (v1.getY() + v2.getY())/2;
    double coordZ = (v1.getZ() + v2.getZ())/2;

    Vertex3D v =  Vertex3D(coordX, coordY, coordZ);
    return v;
}

/**
 * @brief Segmenter::centroidDistance
 * @param f1 index of the first triangle
 * @param f2 index of the second triangle
 * @return approximation of the geodesic distance between
 *      the two triangle barycenters
 */
float Segmenter::centroidDistance(int f1, int f2){

    //first triangle
    Triangle t1 = mesh.GetTopSimplex(f1);

    int indexT, indV1;
    Vertex3D v1, v2, halfP;
    indexT = -1;

    for(int ii=0;ii<3;ii++){

        if(t1.TT(ii)==f2){
            indexT=f2;
            break;
        }
    }
    //There's something wrong, f1 and f2 are not adjacent
    if(indexT < 0){
        cout<<"Not valid faces"<<endl;
        exit(-1);
    }

    // Common edge and vertices
    Edge* shared = t1.TE((indexT)%3);
    indV1 = shared->EV(0);
    v1 = mesh.GetVertex(indV1);
    indV1 = shared->EV(1);
    v2 = mesh.GetVertex(indV1);

    halfP = halfPoint(v1, v2);

    //Corresponding centroids
    Vertex3D C1 = facesCentroids.at(f1);
    Vertex3D C2 = facesCentroids.at(f2);

    return C1.distance(halfP) + halfP.distance(C2);
}

/**
 * @brief Segmenter::buildFaceDistances
 * @return a hash map with the approximate geodesic distances between
 *      pairs of adjacent triangles
 */
std::unordered_map<edgekey, float> Segmenter::buildFaceDistances(){

    std::unordered_map<edgekey, float> FD;

    for(unsigned int ii=0; ii<mesh.GetNumberOfTopSimplexes(); ii++){

        Triangle T = mesh.GetTopSimplex(ii);

        for(int jj=0; jj<3; jj++){

            int f2 = T.TT(jj);
            edgekey ek = getKey(ii, f2);

            //If it is not already in the structure
            if(!FD.count(ek) && f2 >= 0){

                float dist = centroidDistance(ii, f2);
                assert(dist > 0);
                FD[ek] = dist;
            }
        }
    }
    return FD;
}

/**
 * @brief Segmenter::buildAngleDistances
 * @return a hash map with the angular distance between
 *      pairs of adjacent triangles
 */
std::unordered_map<edgekey, float> Segmenter::buildAngleDistances(){

    std::unordered_map<edgekey, float> aDistances;

    for(int ii=0;ii<mesh.GetNumberOfTopSimplexes();ii++){
        Triangle T=mesh.GetTopSimplex(ii);

        for(int jj=0;jj<3;jj++){
            int f2=T.TT(jj);
            edgekey ek=getKey(ii,f2);
            if(aDistances.count(ek)==0 && f2>=0){
                float diffCenters[3];

                //Vertices of the shared edge
                faceind vi1=T.TV((jj+1)%3);
                faceind vi2=T.TV((jj+2)%3);

                Vertex3D v=mesh.GetVertex(vi1);
                Vertex3D w=mesh.GetVertex(vi2);

                double diffx=v.getX()-w.getX();
                double diffy=v.getY()-w.getY();
                double diffz=v.getZ()-w.getZ();

                //Length of shared edge
                double balanceF=sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

                float normF[3];
                Normals N=norms.at(ii);
                float eta;

                //Normal vector for face ii
                normF[0]=N.getNx();
                normF[1]=N.getNy();
                normF[2]=N.getNz();

                Vertex3D C1=facesCentroids.at(ii);
                Vertex3D C2=facesCentroids.at(f2);

                diffCenters[0]=(C2.getX()-C1.getX());
                diffCenters[1]=(C2.getY()-C1.getY());
                diffCenters[2]=(C2.getZ()-C1.getZ());

                Normals nn;
                nn.setNX(diffCenters[0]);
                nn.setNY(diffCenters[1]);
                nn.setNZ(diffCenters[2]);

                if(N.DotProd(nn) > 0){ //concave
                    eta=1.0;
                }
                else{ //convex
                    eta=etaconvex;
                }

                float dot=N.DotProd(norms.at(f2));
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
 * @brief Segmenter::buildGlobalDistances builds the hash map which stores the combined
 *          distance measure between pairs of adjacent triangles
 */
void Segmenter::buildGlobalDistances(){

    for(unsigned int ii=0; ii<mesh.GetNumberOfTopSimplexes(); ii++){

        Triangle T = mesh.GetTopSimplex(ii);

        for(int jj=0; jj<3; jj++){

            int f2=T.TT(jj);
            edgekey ek = getKey(ii, f2);

            if(f2 >= 0 && !globalDistances.count(ek)){

                float dist = (faceDistances[ek] + alpha*angleDistances[ek])/BBDiagonal;
                globalDistances[ek] = dist;
            }
        }
    }
}

/**
 * @brief Segmenter::expansionStep performs a Dijkstra-based expansion
 */
void Segmenter::expansionStep(){

    // initialize graph distance to infinity
    for(unsigned int ii=0;ii<facesCentroids.size();ii++){
        outputDijkstra[ii]=FLT_MAX;
    }

    // initialize clusterIndex expanded count
    for(int i=0;i<facesCentroids.size();i++){
        clusterIndex[i]=-1; // -1 indicates unassigned face
        if(debugMode)
            expanded[i]=0;
    }

    //Time of a step (initialize)
    double accTime=0.0;

    Timer TM2;

    for(unsigned int ii=0;ii<regionCentroids.size();ii++) //initializa the region centroids
        outputDijkstra[regionCentroids[ii]]=0.0;    //  distance to centeroid is always 0


    for(unsigned int ii=0;ii<regionCentroids.size();ii++){
        TM2.start();
        //Visit the triangles with a Dijkstra-based algorithm starting from a region centroid
        expandSeed(regionCentroids[ii], ii);
        if(debugMode)
            expanded[regionCentroids[ii]]++;
        TM2.stop(); //stop the timer
        accTime += TM2.getElapsedTimeInMilliSec();
    }
    double avg=accTime/regionCentroids.size();
    runningTime += avg;

    //Calculate the minimum, maximum and average number of expansions (DEBUG PURPOSES)
    if(debugMode){
        int auxMin=mesh.GetNumberOfTopSimplexes();
        int auxMax=-1;
        float AvgExpansions=0.0;
        int ctrzero=0;
        for(int i=0;i<mesh.GetNumberOfTopSimplexes();i++){
            AvgExpansions+=expanded[i];
            if(expanded[i]<auxMin)
                auxMin=expanded[i];
            if(expanded[i]>auxMax)
                auxMax=expanded[i];
            if(expanded[i]==0){
                ctrzero++;
            }
        }
        AvgExpansions /= mesh.GetNumberOfTopSimplexes();

        cout<<"Min visits = "<<auxMin<<" Max visits = "<<auxMax<<" average = "<<AvgExpansions<<endl;
        cout<<"There are "<<ctrzero<<" unvisited faces"<<endl;
    }

    //Repair step to re-assign faces whose index is -1
    while(!CheckClusterIndex()){
        int violator;
        for(int vv=0;vv<facesCentroids.size();vv++){
            if(/*expanded[vv]==0 || */clusterIndex[vv]<0){
                violator=vv;
                break;
            }
        }
        if(debugMode)
            cout<<"Violator "<<violator<<endl;
        if(violator < facesCentroids.size())
            expandSeed(violator, NCluster++);
    }
    assert(CheckClusterIndex());
}

/**
 * @brief Segmenter::expandSeed expands a single centroid
 * @param indexT index of the triangle which is the center of the current region
 * @param newind index of the current region
 */
void Segmenter::expandSeed(int indexT, int newind){

    priority_queue<pointDist, vector<pointDist>, compare> Q;
    std::unordered_set<faceind> visited;

    if(debugMode){
        if(expanded[indexT]==0)
            expanded[indexT]++;
    }

    pointDist seed;
    seed.indexP=indexT;
    seed.distanceP=0.0;
    Q.push(seed);

    Vertex3D rc = facesCentroids.at(indexT);
    clusterIndex[indexT]=newind;

    regionCentroids[newind]=indexT;

    while(!Q.empty()){

        pointDist actual=Q.top();
        Q.pop();

        if(rc.distance(facesCentroids.at(actual.indexP))/BBDiagonal <= maxD*timesR){
            int neigh;
            Triangle T = mesh.GetTopSimplex(actual.indexP);

            for(int jj=0;jj<3;jj++){
                neigh=T.TT(jj);
                if(visited.count(neigh)==0 && neigh>=0){

                    float newdist=actual.distanceP + globalDistances[getKey(actual.indexP, neigh)];
                    if(debugMode)
                        expanded[neigh]++;
                    if(newdist < outputDijkstra[neigh]){
                        outputDijkstra[neigh]=newdist;
                        clusterIndex[neigh]=newind;

                        pointDist pd;
                        pd.indexP=neigh;
                        pd.distanceP=newdist;
                        Q.push(pd);
                    }
                }
            }
            visited.insert(actual.indexP);
        }
    }
}

/**
 * @brief Segmenter::placeSeeds Initialization if we pass the number of desired region as parameter
 * @param index index of the first centroid
 */
void Segmenter::placeSeeds(int index){

    regionCentroids[0] = index;
    clusterIndex[index]=0;

    double *distFromSeeds = new double[mesh.GetNumberOfTopSimplexes()];

    int count=1;
    cout<<"Initializing..."<<endl;
    while(regionCentroids.size()<NCluster){
        for(int iterFaces = 0; iterFaces < mesh.GetNumberOfTopSimplexes(); iterFaces++){
            distFromSeeds[iterFaces] = 0.0;
            if(regionCentroids.count(iterFaces)==0){

                distFromSeeds[iterFaces] = facesCentroids.at(iterFaces).distance(facesCentroids.at(regionCentroids[0]));
                for(int kk=1;kk<regionCentroids.size();kk++){

                    if(distFromSeeds[iterFaces] > facesCentroids.at(iterFaces).distance(facesCentroids.at(regionCentroids[kk])))
                        distFromSeeds[iterFaces] = facesCentroids.at(iterFaces).distance(facesCentroids.at(regionCentroids[kk]));

                }
            }
        }
        faceind indOfM = indexOfMax(distFromSeeds);
        clusterIndex[indOfM] = count;
        regionCentroids[count++] = indOfM;
    }
    if(debugMode)
        cout<<"Done: regions "<<regionCentroids.size()<<endl;
}

/**
 * @brief Segmenter::initializationGrid initialize the segmentation dividing the object with a regular grid
 */
void Segmenter::initializationGrid(){

    std::unordered_map<gridkey, int> hmap; //for already assigned maps
    int id=0; //index of current region
    float denominator = 2*maxD*BBDiagonal;

    for(unsigned int ii=0;ii<facesCentroids.size();ii++){

        Vertex3D actualF = facesCentroids.at(ii);
        //find which region it belongs to
        int i = floor(actualF.getX()/denominator);
        int j = floor(actualF.getY()/denominator);
        int k = floor(actualF.getZ()/denominator);

        //if it is the first time we encounter the region
        if(hmap.count(getGridKey(i,j,k))==0){
            regionCentroids[id]=ii;
            hmap[getGridKey(i,j,k)] = id++;
        }

        //assign a region index to the triangle
        clusterIndex[ii]=hmap[getGridKey(i,j,k)];
    }
}

/**
 * @brief Segmenter::floodInitialization initialization with a Dijkstra-based expansion
 * @param indexT index of the startin point
 */
void Segmenter::floodInitialization(int indexT){

    //Index of the actual region
    int actualInd=0;

    //Array to keep track of expanded faces (DEBUG PURPOSES)
    if(debugMode){

        int *expanded = new int[mesh.GetNumberOfTopSimplexes()];
        for(int i=0;i<mesh.GetNumberOfTopSimplexes();i++){
            expanded[i]=0;
        }
        expanded[indexT]++; //Mark the seed as expanded
    }

    // Initialize the distance graph to infinity for each triangle
    for(unsigned int ii=0;ii<mesh.GetNumberOfTopSimplexes();ii++){
        outputDijkstra[ii]=FLT_MAX;
    }
    outputDijkstra[indexT]=0; //Distance of the seed is 0

    // Initialize index of each face to -1 (unassigned)
    for(unsigned int ii=0;ii<mesh.GetNumberOfTopSimplexes();ii++)
        clusterIndex[ii]=-1;

    //Queue that stores all the triangles ordered with respect to their Euclidean distance from the seed
    //we had it to avoid having unassigned faces at the end of the initialization
    cout<<"Enter in the first loop"<<endl;
    for(unsigned int a=0;a<mesh.GetNumberOfTopSimplexes();a++){
        pointDist pp;
        pp.indexP=a;
        if(a==indexT)
            pp.distanceP=0.0;
        else
            pp.distanceP=facesCentroids.at(a).distance(facesCentroids.at(indexT));

        globalQ.push(pp);
    }

    // Main loop
    while(!globalQ.empty()){

        // Remove faces already assigned to some cluster
        while(clusterIndex[globalQ.top().indexP]>=0 && globalQ.size()>1)
            globalQ.pop();

        //Take the first element of the queue
        pointDist center=globalQ.top();
        globalQ.pop();

        if(center.distanceP==FLT_MAX){ //Unconnected or unvisited face
            cout<<"Unconnected or unvisited"<<endl;
            break;
        }

        if(debugMode){
            if(expanded[regionCentroids[actualInd]]==0)
                expanded[regionCentroids[actualInd]]++;
        }
        expandSeed(center.indexP, actualInd);
        actualInd++;
    }
}

/**
 * @brief Segmenter::getBBDiagonal sets lower left and upper right corners of the mesh bounding box
 *          and its diagonal length
 */
void Segmenter::getBBDiagonal(){

    //Minimum
    minCoords.setX(mesh.GetVertex(0).getX());
    minCoords.setY(mesh.GetVertex(0).getY());
    minCoords.setZ(mesh.GetVertex(0).getZ());

    //Maximum
    maxCoords.setX(mesh.GetVertex(0).getX());
    maxCoords.setY(mesh.GetVertex(0).getY());
    maxCoords.setZ(mesh.GetVertex(0).getZ());

    for(int ii=0; ii<mesh.GetNumVertices(); ii++){

        Vertex3D vv=mesh.GetVertex(ii);

        //Update if necessary
        if(vv.getX()<minCoords.getX())
            minCoords.setX(vv.getX());
        if(vv.getY()<minCoords.getY())
            minCoords.setY(vv.getY());
        if(vv.getZ()<minCoords.getZ())
            minCoords.setZ(vv.getZ());

        if(vv.getX()>maxCoords.getX())
            maxCoords.setX(vv.getX());
        if(vv.getY()>maxCoords.getY())
            maxCoords.setY(vv.getY());
        if(vv.getZ()>maxCoords.getZ())
            maxCoords.setZ(vv.getZ());
    }

    //diagonal length
    BBDiagonal = maxCoords.distance(minCoords);
}

/**
 * @brief Segmenter::CheckClusterIndex
 * @return true if every triangle is assigned to some cluster (its index is >=0)
 */
bool Segmenter::CheckClusterIndex(){
    bool ret=true;

    for(int a=0;a<mesh.GetNumberOfTopSimplexes();a++){
        if(clusterIndex[a] < 0){
            ret=false;
        }
    }
    return ret;
}

/**
 * @brief minArrayIndex
 * @param array
 * @return the index of the minimum value into the array
 */
int minArrayIndex(float array[3]){
    if(array[0]<array[1] && array[0]<array[2])
        return 0;
    else if(array[1]<array[2])
        return 1;
    else
        return 2;
}

/**
 * @brief Segmenter::updateCenters takes centroids as close as possible to the center of the region
 * @return true if no centroid has moved from the previous iteration
 */
bool Segmenter::updateCenters(){

    // differences in each region
    int differences[NCluster];

    for(int k=0;k<NCluster;k++)
        differences[k]=0;

    bool any_moves=false;

    // stores the old centroids
    std::unordered_map<edgekey, int> olds;
    if(debugMode)
        cout<<"Update"<<endl;

    for(unsigned int ii=0;ii<regionCentroids.size();ii++)
        olds[ii]=regionCentroids[ii];

    // Iterate on each one of the regions
    for(unsigned int ii=0;ii<regionCentroids.size();ii++){

        //Area of a region
        double regionArea=0.0;
        for(unsigned int ff=0;ff<facesCentroids.size();ff++){
            if(clusterIndex[ff]==ii)
                regionArea += mesh.TriangleArea(ff);
        }
        double xc=0.0;
        double yc=0.0;
        double zc=0.0;

        // Area-weighted average
        for(unsigned int ff=0;ff<facesCentroids.size();ff++){
            if(clusterIndex[ff]==ii){
                Vertex3D current=facesCentroids.at(ff);
                double ta=mesh.TriangleArea(ff);
                xc += ta*current.getX();
                yc += ta*current.getY();
                zc += ta*current.getZ();
            }
        }
        //Average of the coordinates of the baricenters of the region
        xc /= regionArea;
        yc /= regionArea;
        zc /= regionArea;

        Vertex3D nc;
        nc.setX(xc);
        nc.setY(yc);
        nc.setZ(zc);

        double actualD=nc.distance(facesCentroids.at(regionCentroids[ii]));

        //Centroid closest to nc
        for(unsigned int ff=0;ff<facesCentroids.size();ff++){
            if(clusterIndex[ff]==ii){
                if(nc.distance(facesCentroids.at(ff)) < actualD && ff != regionCentroids[ii]){
                    actualD=nc.distance(facesCentroids.at(ff));
                    any_moves=true;
                    regionCentroids[ii]=ff;
                    differences[ii]++;
                }
            }
        }
    }

    // count the differences
    int ctrdiff=0;
    int totaldiff=0;
    for(int k=0;k<NCluster;k++){
        if(differences[k])
            ctrdiff++;
        totaldiff+=differences[k];
    }
    if(debugMode){
        cout<<"Regions "<<NCluster<<endl;
        cout<<"Differences in "<<ctrdiff<<" regions, total = "<<totaldiff<<endl;
    }
    return any_moves;
}

/**
 * @brief Segmenter::writeSegmOnFile
 * @param fileout the file in which the segmentation is saved
 * @return 0 if everything is ok, -1 otherwise
 */
int Segmenter::writeSegmOnFile(string fileout){

    //open file in write mode
        ofstream ofs;
        ofs.open(fileout.c_str());

        time_t now;
        time(&now);

        struct tm * localT=localtime(&now);

        if(ofs.is_open()){

            //if we want an header (visualization purposes)
            if(putHeader){
                ofs << "#";
                ofs << asctime(localT) << endl;
                ofs << "#mesh input=" << filename <<endl;
                ofs << "#radius=" << maxD <<endl;
                ofs << "#alpha=" << alpha <<endl;
                ofs << "#eta convex=" <<etaconvex <<endl;
                //number of iterations
                ofs << "#iterations: "<< iters << endl;
                //running time of the segmentation
                ofs << "#time: "<< millisecs << endl << endl;
            }
            for(int a=0;a<facesCentroids.size();a++){
                ofs << clusterIndex[a] << endl;
            }

            ofs.close();

            return 0;
        }
        else{
            cout<<"Unable to write on file"<<endl;
            return -1;
        }
}
