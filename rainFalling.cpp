#include "segmenter.h"

vector<int> Segmenter::getMinAndCLuster(vector<Cluster *> *clusters){

    vector<bool> visited = vector<bool>(clusters->size(), false);
    vector<int> minimum_labels = vector<int>(clusters->size(), -1);
    int label=0;

    for(unsigned int ii=0; ii<clusters->size(); ii++){

        if(visited[ii])
            continue;

        queue<int> plateau;
        plateau.push(ii);
        Cluster* clst = new Cluster();
        clst->vertices.insert(ii);

        while (!plateau.empty()) {
           int vertex = plateau.front();
           visited[vertex]=true;
           plateau.pop();

           vector<int> link = mesh.VV(ii);
           for(int jj=0; jj<link.size(); jj++){

               if(mesh.getVertex(link.at(jj)).getZ() == mesh.getVertex(vertex).getZ() && !visited[link.at(jj)]){
                   clst->vertices.insert(link.at(jj));
                   visited[link.at(jj)]=true;
                   plateau.push(link.at(jj));
               }

               if(mesh.getVertex(link.at(jj)).getZ() < mesh.getVertex(vertex).getZ())
                   clst->tipo=NORMAL;
           }
        }

        if(clst->vertices.size() > 1){
            for(set<int>::iterator it=clst->vertices.begin(); it!=clst->vertices.end(); it++)
                clusters->at(*it)=clst;
        }
        if(clst->tipo==MINIMUM){
            for(set<int>::iterator it=clst->vertices.begin(); it!=clst->vertices.end(); it++)
                minimum_labels[*it]=label;
            label++;
        }
    }
    for(int i=0; i<visited.size(); i++)
        assert(visited[i]);

    return minimum_labels;
}


int Segmenter::find_next(int v, vector<int> *minima, vector<Cluster *> *clusters){

    if(clusters->at(v) != NULL){
        int min=v;

        for(set<int>::iterator it=clusters->at(v)->vertices.begin(); it!=clusters->at(v)->vertices.end(); it++){
            vector<int> link=mesh.VV(*it);
            for(int i=0;i<link.size();i++){
                if(mesh.getVertex(min).getZ() > mesh.getVertex(link[i]).getZ())
                    min=link.at(i);
            }
        }
        assert(min != v && mesh.getVertex(min).getZ() < mesh.getVertex(v).getZ());

        int label=minima->at(min);
        if(minima->at(min)==-1)
            label=find_next(min, minima, clusters);
        assert(label!=-1 && label==minima->at(min));

        for(set<int>::iterator it = clusters->at(v)->vertices.begin(); it!=clusters->at(v)->vertices.end(); it++)
            minima->at(*it)=minima->at(min);
        return minima->at(min);
    }

    else{
        vector<int> link = mesh.VV(v);
        int min=v;

        for(int i=0; i<link.size(); i++){
            if(mesh.getVertex(min).getZ() > mesh.getVertex(link.at(i)).getZ())
                min=link.at(i);
        }
        assert(min != v);

        int label=minima->at(min);
        if(minima->at(min)==-1){
            label=find_next(min, minima, clusters);
        }
        assert(label!=-1 && label==minima->at(min));
        minima->at(v)=minima->at(min);

        return minima->at(min);
    }
}

void Segmenter::descend_paths(vector<int> *minima, vector<Cluster *> *clusters){

    for(int i=0; i<clusters->size(); i++){
        if(clusters->at(i) != NULL && clusters->at(i)->tipo==MINIMUM && minima->at(i)==-1)
            int end_label=this->find_next(i, minima, clusters);
    }

    for(int i=0; i<mesh.getNumVertex(); i++){
        if(minima->at(i) ==-1)
            int end_label=find_next(i, minima, clusters);
    }
}

int* Segmenter::label_triangles(vector<int> *minima){

    //vector<int> label_tri(mesh.getTopSimplexesNum(), 0);
    int* label_tri = new int[mesh.getTopSimplexesNum()];
    for(int ii=0;ii<mesh.getTopSimplexesNum();ii++)
        label_tri[ii]=0;

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        int lower=mesh.getTopSimplex(ii).TV(0);
        for(int jj=0; jj<3; jj++){
            if(mesh.getVertex(mesh.getTopSimplex(ii).TV(jj)).getZ() < mesh.getVertex(lower).getZ())
                lower=mesh.getTopSimplex(ii).TV(jj);
        }
        label_tri[ii]=minima->at(lower);
    }
    return label_tri;
}

void Segmenter::loadWatershed(){

    openMeshFile(this->filename);
    vector<Cluster*> clusters = vector<Cluster*>(mesh.getNumVertex(), NULL);
    this->clusterIndex = new int[mesh.getTopSimplexesNum()];

    vector<int> minimum_labels = getMinAndCLuster(&clusters);
    descend_paths(&minimum_labels, &clusters);

    clusterIndex = label_triangles(&minimum_labels);
    countRegions();
}
