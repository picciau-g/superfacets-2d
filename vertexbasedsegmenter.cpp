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
}

void VertexBasedSegmenter::loadMesh(){

    openMeshFile(this->filename);
    this->centerMesh = centerCoordinate();

    faceAreas = new float[mesh.getTopSimplexesNum()];
    clusterIndex = new int[mesh.getNumVertex()];

    for(unsigned int ii=0; ii<mesh.getTopSimplexesNum(); ii++){
        Triangle T = mesh.getTopSimplex(ii);
        Normals n = Normals(mesh.getVertex(T.TV(0)), mesh.getVertex(T.TV(1)), mesh.getVertex(T.TV(2)));
        norms.push_back(n);
    }
    setAreas();
    getBBDiagonal();

    for(int ii=0; ii<mesh.getNumVertex(); ii++)
        clusterIndex[ii]=-1;

    vertexDistances = buildVertexDistances();
    functionVDistances = buildFunctionVDistances();
    buildGlobalDistances();

    vertexDistances.erase(vertexDistances.begin(), vertexDistances.end());
    functionVDistances.erase(functionVDistances.begin(), functionVDistances.end());

    startSeg();
}

void VertexBasedSegmenter::startSeg(){

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


void VertexBasedSegmenter::Segmentation(){

    while(updateCenters())
        expansionStep();

    cout<<"Ok, converged"<<endl;
}

void VertexBasedSegmenter::expansionStep(){


    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){
        outputDijkstra[ii]=FLT_MAX;
        clusterIndex[ii]=-1;
    }

    for(unsigned int ii=0; ii<regionCentroids.size(); ii++)
        outputDijkstra[regionCentroids[ii]]=0.0;    //  distance to centeroid is always 0

    for(unsigned int ii=0; ii<regionCentroids.size(); ii++){
        expandSeed(regionCentroids[ii], ii);
    }

    while (!CheckClusterIndex()) {
        int violator;
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
    std::tr1::unordered_set<vertexind> visited;

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

        if(rc.distance(mesh.getVertex(actual.indexP))/BBDiagonal <= maxD){
            int neigh;
            Vertex3D V = mesh.getVertex(actual.indexP);

            vector<int> VV = mesh.VV(actual.indexP);
            for(unsigned int ii=0; ii<VV.size(); ii++){
                neigh = VV.at(ii);

                if(neigh >= 0 && !visited.count(neigh)){

                    float newdist = actual.distanceP + globalDistances[getKey(actual.indexP, neigh)];
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
    std::tr1::unordered_map<edgekey, int> olds;

    for(unsigned int ii=0;ii<regionCentroids.size();ii++){

        double xc=0.0, yc=0.0, zc=0.0;
        int cardinality=0;
        for(unsigned int jj=0; ii<mesh.getNumVertex(); jj++){

            if(clusterIndex[jj]==ii){
                Vertex3D V = mesh.getVertex(jj);
                xc += V.getX();
                yc += V.getY();
                zc += V.getZ();

                cardinality++;
            }
        }
        xc /= cardinality;
        yc /= cardinality;
        zc /= cardinality;

        Vertex3D nc(xc, yc, zc);

        double actualD = nc.distance(mesh.getVertex(regionCentroids[ii]));

        for(unsigned int jj=0; jj<mesh.getNumVertex(); jj++){

            if(clusterIndex[jj]==ii){

                if(nc.distance(mesh.getVertex(jj)) < actualD && jj != regionCentroids[jj]){
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

std::tr1::unordered_map<edgekey, float> VertexBasedSegmenter::buildVertexDistances(){

    std::tr1::unordered_map<edgekey, float> VD;

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D V = mesh.getVertex(ii);
        vector<int> VNeighbors = mesh.VV(ii);

        for(unsigned int jj=0; jj<VNeighbors.size(); jj++){

            int n = VNeighbors.at(jj);
            edgekey ek = getKey(ii, n);

            if(!VD.count(ek) && n>=0){

                float dist = V.distance(mesh.getVertex(n));
                assert(dist>0);
                VD[ek] = dist;
            }
        }
    }
    return VD;
}

std::tr1::unordered_map<edgekey, float> VertexBasedSegmenter::buildFunctionVDistances(){

    std::tr1::unordered_map<edgekey, float> fdDistances;

    for(int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D V = mesh.getVertex(ii);
        vector<int> VNeighbors = mesh.VV(ii);

        for(unsigned int jj=0; jj<VNeighbors.size(); jj++){

            int nv = VNeighbors.at(jj);
            edgekey ek = getKey(ii, nv);

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
            edgekey ek = getKey(ii, n);

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

    FILE* f = fopen(curvFile.c_str(), 'r');

    int num;
    float fv;
    fscanf(f, "%d", &num);

    for(int ii=0; ii<num; ii++){
        fscanf(f, "%f", &fv);
        functionValue.push_back(fv);
    }

    cout<<"read "<<functionValue.size()<<"elements"<<endl;
}
