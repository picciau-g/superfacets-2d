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

}

void VertexBasedSegmenter::expandSeed(int a, int b){

}

bool VertexBasedSegmenter::updateCenters(){

}

std::tr1::unordered_map<edgekey, float> VertexBasedSegmenter::buildVertexDistances(){

}

std::tr1::unordered_map<edgekey, float> VertexBasedSegmenter::buildFunctionVDistances(){

}

void VertexBasedSegmenter::buildGlobalDistances(){

}

void VertexBasedSegmenter::getBBDiagonal(){

}
