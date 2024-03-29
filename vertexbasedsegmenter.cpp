#include "vertexbasedsegmenter.h"

VertexBasedSegmenter::VertexBasedSegmenter()
{

}

void VertexBasedSegmenter::openMeshFile(string mname){

    mesh = Mesh<Vertex3D, Triangle>();
    QString qmn = QString::fromStdString(mname);
    QStringList qsl = qmn.split(".");

    if(!qsl.back().compare("tri"))
        Reader::readMeshFile(mesh, mname);
    else if(!qsl.back().compare("off"))
        Reader::readOFFMesh(mesh, mname);

    else{
        cout<<"Not a valid file format (it must be .tri or .off)"<<endl;
        exit(0);
    }
    mesh.build();

    cout<<"Build mesh with "<<mesh.getNumVertex()<<" vertices and "<<mesh.getTopSimplexesNum()<<" triangles "<<endl;
}

void VertexBasedSegmenter::loadMesh(){

    openMeshFile(this->filename);
    cout<<"Loading "<<this->filename<<endl;
    this->centerMesh = centerCoordinate();

    faceAreas = new float[mesh.getTopSimplexesNum()];
    clusterIndex = new int[mesh.getNumVertex()];
    cout<<"found center"<<endl;

    for(unsigned int ii=0; ii<mesh.getTopSimplexesNum(); ii++){
        Triangle T = mesh.getTopSimplex(ii);
        Normals n = Normals(mesh.getVertex(T.TV(0)), mesh.getVertex(T.TV(1)), mesh.getVertex(T.TV(2)));
        norms.push_back(n);
    }
    cout<<"set normals"<<endl;
    setAreas();
    getBBDiagonal();
    cout<<"Diag "<<this->BBDiagonal<<endl;

    double auxRad = sqrt(mesh.MArea()/(NCluster*M_PI));
    this->maxD = auxRad/BBDiagonal;

    openCurvatureFile(fieldfilename);
    cout<<"Loaded function"<<endl;

    for(int ii=0; ii<mesh.getNumVertex(); ii++)
        clusterIndex[ii]=-1;

    cout<<"Before vertices"<<endl;
    vertexDistances = buildVertexDistances();
    cout<<"Vrtices built"<<endl;
    functionVDistances = buildFunctionVDistances();
    cout<<"function built"<<endl;
    buildGlobalDistances();
    cout<<"global built"<<endl;

    vertexDistances.erase(vertexDistances.begin(), vertexDistances.end());
    functionVDistances.erase(functionVDistances.begin(), functionVDistances.end());

    cout<<"All built, "<<mesh.getNumVertex()<<" vertices and "<<mesh.getTopSimplexesNum()<<" triangles"<<endl;

    //startSeg();
}

Vertex3D VertexBasedSegmenter::centerCoordinate(){

    Vertex3D res;
    float sumX=0.0; float sumY=0.0; float sumZ=0.0;

    for(int ii=0;ii<mesh.getNumVertex();ii++){

        sumX += mesh.getVertex(ii).getX();
        sumY += mesh.getVertex(ii).getY();
        sumZ += mesh.getVertex(ii).getZ();
    }
    res.setX(sumX/mesh.getNumVertex());
    res.setY(sumY/mesh.getNumVertex());
    res.setZ(sumZ/mesh.getNumVertex());

    return res;
}

int VertexBasedSegmenter::nearestVertex(){

    float distance = centerMesh.distance(mesh.getVertex(0));
    int nearestSoFar = 0;

    for(int v=0; v<mesh.getNumVertex(); v++){

        if(distance>centerMesh.distance(mesh.getVertex(v))){
            distance = centerMesh.distance(mesh.getVertex(v));
            nearestSoFar = v;
        }
    }
    return nearestSoFar;
}

Vertex3D VertexBasedSegmenter::halfPoint(Vertex3D v1, Vertex3D v2){

    double coordX = (v1.getX() + v2.getX())/2;
    double coordY = (v1.getY() + v2.getY())/2;
    double coordZ = (v1.getZ() + v2.getZ())/2;

    Vertex3D v =  Vertex3D(coordX, coordY, coordZ);
    return v;
}

void VertexBasedSegmenter::placeSeeds(int index){

    regionCentroids[0] = index;
    clusterIndex[index] = 0;
    int count = 1;

    double *distFromSeeds = new double[mesh.getNumVertex()];

    while (regionCentroids.size() < NCluster) {

        for(int iterVerts = 0; iterVerts<mesh.getNumVertex(); iterVerts++){
            distFromSeeds[iterVerts] = 0.0;

            if(!regionCentroids.count(iterVerts)){

                distFromSeeds[iterVerts] = mesh.getVertex(iterVerts).distance(mesh.getVertex(regionCentroids[0]));
                for(int k=1; k<regionCentroids.size(); k++){
                    if(distFromSeeds[iterVerts] > mesh.getVertex(iterVerts).distance(mesh.getVertex(regionCentroids[k])))
                        distFromSeeds[iterVerts] = mesh.getVertex(iterVerts).distance(mesh.getVertex(regionCentroids[k]));
                }

            }
        }

        vertexind furthest = indexOfMax(distFromSeeds);
        clusterIndex[furthest] = count;
        regionCentroids[count++] = furthest;
        //cout<<"Reg "<<count-1<<" vertex "<<regionCentroids[count-1]<<endl;
    }
    cout<<"Placed "<<regionCentroids.size()<<" centroids"<<endl;
}


void VertexBasedSegmenter::Segmentation(){

    placeSeeds(0);
    expansionStep();
    cout<<"Reg "<<regionCentroids.size()<<endl;
    int iterNumber = 0;

    cout<<"Initialized, now start loop..."<<endl;
    while(updateCenters() && iterNumber<20){
        expansionStep();
        cout<<"Done iteration "<<iterNumber++<<"..."<<endl;
    }

    cout<<"Ok, converged (regions "<<regionCentroids.size()<<")"<<endl;
}

void VertexBasedSegmenter::expansionStep(){


    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){
        outputDijkstra[ii]=FLT_MAX;
        clusterIndex[ii]=-1;
    }

    //cout<<"Max "<<maxD<<endl;

    for(unsigned int ii=0; ii<NCluster; ii++)
        outputDijkstra[regionCentroids[ii]]=0.0;    //  distance to centeroid is always 0

    for(unsigned int ii=0; ii<NCluster; ii++){
        expandSeed(regionCentroids[ii], ii);
    }

    while (!CheckClusterIndex()) {
        int violator;
        cout<<"Violator"<<endl;
        for(int vv=0; vv<mesh.getNumVertex(); vv++){
            if(clusterIndex[vv]<0){
                violator=vv;
                break;
            }
        }
        if(violator<mesh.getNumVertex())
            expandSeed(violator, NCluster++);
    }
}

void VertexBasedSegmenter::expandSeed(int indexV, int newInd){

    priority_queue<pointDist, vector<pointDist>, compare> Q;
    std::unordered_set<vertexind> visited;

    pointDist seed;
    seed.indexP = indexV;
    seed.distanceP = 0.0;
    Q.push(seed);

    Vertex3D rc = mesh.getVertex(indexV);
    clusterIndex[indexV] = newInd;

    while(!Q.empty()){

        pointDist actual = Q.top();
        Q.pop();

        if(actual.indexP > mesh.getNumVertex()){
            cout<<"Index exceed"<<endl;
            continue;
        }

        if(rc.distance(mesh.getVertex(actual.indexP))/BBDiagonal <= 2*maxD){
            int neigh;
            Vertex3D V = mesh.getVertex(actual.indexP);

            vector<int> VV = mesh.VV(actual.indexP);
            for(unsigned int ii=0; ii<VV.size(); ii++){
                neigh = VV.at(ii);

                if(neigh >= 0 && !visited.count(neigh)){

                    float newdist = actual.distanceP + globalDistances[getVertexKey(actual.indexP, neigh)];
                    if(newdist < outputDijkstra[neigh]){
                        outputDijkstra[neigh] = newdist;
                        clusterIndex[neigh] = newInd;

                        pointDist pd;
                        pd.indexP = neigh;
                        pd.distanceP = newdist;
                        Q.push(pd);
                    }
                }
            }
            visited.insert(actual.indexP);
        }
    }
}

bool VertexBasedSegmenter::updateCenters(){

    int differences[NCluster];

    for(unsigned int ii=0; ii<NCluster; ii++)
        differences[ii]=0;

    bool any_moves=false;
    std::unordered_map<edgekey, int> olds;

    for(unsigned int ii=0;ii<NCluster;ii++){

        double xc=0.0, yc=0.0, zc=0.0;
        int cardinality=0;
        //cout<<"Iter "<<ii<<endl;
        for(unsigned int jj=0; jj<mesh.getNumVertex(); jj++){

            if(clusterIndex[jj]==ii){
                Vertex3D V = mesh.getVertex(jj);
                xc += V.getX();
                yc += V.getY();
                zc += V.getZ();

                cardinality++;
            }
        }

        if(cardinality>0){
            xc /= cardinality;
            yc /= cardinality;
            zc /= cardinality;
        }

        Vertex3D nc(xc, yc, zc);

        double actualD = nc.distance(mesh.getVertex(regionCentroids[ii]));

        for(unsigned int jj=0; jj<mesh.getNumVertex(); jj++){

            if(clusterIndex[jj]==ii){

                if(nc.distance(mesh.getVertex(jj)) < actualD && jj != regionCentroids[ii]){
                    actualD = nc.distance(mesh.getVertex(jj));
                    any_moves=true;
                    regionCentroids[ii]=jj;
                    differences[ii]++;
                }
            }
        }
    }

    //count the differences
    int totaldiff=0;
    for(int k=0; k<NCluster; k++){
        totaldiff+=differences[k];
    }
    return any_moves;
}

std::unordered_map<edgekey, float> VertexBasedSegmenter::buildVertexDistances(){

    std::unordered_map<edgekey, float> VD;

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        //cout<<"Loop "<<ii<<endl;
        Vertex3D V = mesh.getVertex(ii);
        //cout<<"V "<<V.getX()<<" "<<V.getY()<<" "<<V.getZ()<<endl;
        vector<int> VNeighbors = mesh.VV(ii);
        //cout<<"Loop "<<ii<<" size "<<VNeighbors.size()<<endl;
        //cout<<"I "<<ii<<endl;

//        for(int jj=0; jj<VNeighbors.size(); jj++)
//            cout<<"J "<<jj<<" N "<<VNeighbors.at(jj)<<" ";
//        cout<<endl;

        if(VNeighbors.size() == 0)
            continue;

        for(unsigned int jj=0; jj<VNeighbors.size(); jj++){

            int n = VNeighbors.at(jj);
            if(n<0 || n>mesh.getNumVertex())
                continue;
            edgekey ek = getVertexKey(ii, n);

            if(!VD.count(ek) && n>=0){

                float dist = V.distance(mesh.getVertex(n));
                //assert(dist>0);
                VD[ek] = dist;
            }
        }
        //cout<<"out of loop "<<ii<<endl;
    }
    cout<<"Returning"<<endl;
    return VD;
}

std::unordered_map<edgekey, float> VertexBasedSegmenter::buildFunctionVDistances(){

    std::unordered_map<edgekey, float> fdDistances;

    for(int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D V = mesh.getVertex(ii);
        vector<int> VNeighbors = mesh.VV(ii);

        for(unsigned int jj=0; jj<VNeighbors.size(); jj++){

            int nv = VNeighbors.at(jj);
            edgekey ek = getVertexKey(ii, nv);

            if(!fdDistances.count(ek) && nv>=0){
                float D = fabs(functionValue.at(ii)-functionValue.at(jj));
                fdDistances[ek] = D;
            }
        }
    }
    return fdDistances;
}

void VertexBasedSegmenter::buildGlobalDistances(){

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        vector<int> VNeighbors = mesh.VV(ii);

        for(unsigned int jj=0; jj<VNeighbors.size(); jj++){

            int n = VNeighbors.at(jj);
            edgekey ek = getVertexKey(ii, n);

            if(n>=0 && !globalDistances.count(ek)){
                float dist = vertexDistances[ek]/BBDiagonal + alpha*functionVDistances[ek];
                globalDistances[ek] = dist;
            }
        }
    }
}

void VertexBasedSegmenter::getBBDiagonal(){

    minCoords.setX(mesh.getVertex(0).getX());
    minCoords.setY(mesh.getVertex(0).getY());
    minCoords.setZ(mesh.getVertex(0).getZ());

    maxCoords.setX(mesh.getVertex(0).getX());
    maxCoords.setY(mesh.getVertex(0).getY());
    maxCoords.setZ(mesh.getVertex(0).getZ());

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D vv = mesh.getVertex(ii);
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

    BBDiagonal = maxCoords.distance(minCoords);
}

int VertexBasedSegmenter::writeSegmOnFile(string filename){

    ofstream ofs;
    ofs.open(filename.c_str());

    if(ofs.is_open()){

        for(int a=0; a<mesh.getNumVertex(); a++)
            ofs << clusterIndex[a] <<endl;

        ofs.close();
        return 0;
    }

    cout<<"Unable to write on file "<<endl;
    return -1;
}


void VertexBasedSegmenter::openCurvatureFile(string curvFile){

    FILE* f = fopen(curvFile.c_str(), "r");

    int num;
    float fv;
    fscanf(f, "%d", &num);

    for(int ii=0; ii<num; ii++){
        fscanf(f, "%f", &fv);
        functionValue.push_back(fv);
    }

    cout<<"read "<<functionValue.size()<<"elements"<<endl;
}


edgekey getVertexKey(vertexind a, vertexind b){

    if(a<=b)
        return (edgekey(a) << 32 | edgekey(b));
    else
        return (edgekey(b) << 32 | edgekey(a));
}
