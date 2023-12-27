#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <assert.h>

#include "Sorting.h"
#include "Edge.h"

/*
True iff v is an endpoint of edge v1-v2
*/
#define is_endpoint(v,v1,v2) ( ( (v)==(v1) ) || ( (v)==(v2) ) )
enum versus { CW=0 , CCW=1 };
using namespace std;

/// A class representing a generic mesh parametrized by the type of top simplexes
template<class V, class T>
class Mesh
{
public:
    ///A constructor method
    Mesh();
    ///A constructor method
    Mesh(const Mesh& orig);
    ///A destructor method
    virtual ~Mesh();
    ///A public method that returns the vertex at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a Vertex&, the vertex at the id-th position in the list
     */
    V& getVertex(int id);
    ///A public method that returns the number of mesh vertices
    /*!
     * \return an integer, representing the number of vertices
     */
    int getNumVertex();
    ///A public method that adds a vertex to the vertices list
    /*!
     * \param v a Vertex& argument, representing the vertex to add
     */
    void addVertex(V& v);
    ///A public method that returns the tetrahedron at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a T&, the tetrahedron at the id-th position in the list
     */
    T& getTopSimplex(int id);
    ///A public method that returns the number of top simplexes
    /*!
     * \return an integer, representing the number of top simplexes
     */
    int getTopSimplexesNum();
    ///A public method that adds a top simplex to the top simplexes list
    /*!
     * \param t a T& argument, representing the top simplex to add
     */
    void addTopSimplex(T& t);
    ///A public method that initializes the space needed by the vertices and top simplexes lists
    /*!
     * \param numV an integer, represents the number of vertices
     * \param numT an integer, represents the number of top simplexes
     */
    void reserveVectorSpace(int numV, int numT);
    void reserveVectorSpace_TopSimplexes(int numT);
    void reserveVectorSpace_Vertices(int numV);

    void build();

    bool half_edge_collapse(int v1, int v2, int t1, int t2);
    bool convex_neighborhood(int v1, int v2, int t1, int t2);


    vector<int> VT(int center);
    vector<Edge*> VE(int center);
    vector<int> VV(int center);

    bool isBoundary(int center);
    bool is_alive(int f){return !removed_triangle[f];}
    bool is_v_alive(int v){return !removed_vertex[v];}
    void reorder_triangulation();


    //compute cosinus of angle formed by 3 vertices
    double cosAngle(int v1, int v2, int v3);
    //compute total area of triangles incident in v
    double FanArea(int v);
    //compute triangle area
    double TArea(int t);
    //computes area of the whole mesh
    double MArea();
    //barycenter
    void Barycenter();

    inline vector<double> retBarycenter(){
        return this->coordsOfB;
    }

    int VIndexInT(int v, int t);
    int NextTAroundV(int t, int v, versus verso);

    //compute the mixed (Voronoi-barycentric) area around vertex v
    double VoronoiBarycentricArea(int v);
    //compute either the Voronoi area (if all angles of t are acute) or the
    //barycentric area (otherwise) of triangle t centered in its vertex v
    double VoronoiBarycentricArea(int v, int t);


    vector<int> ET(Edge e);
    vector<Edge*> EE(Edge e);


    //Utility Function
    template<class C> bool isIntoVector(C c, vector<C> &c_vect);

protected:
    ///A private varible representing the vertices list of the mesh
    vector<V> vertices;
    ///A private varible representing the top simplexes list of the mesh
    vector<T> topSimplexes;

    vector<bool> removed_vertex;
    vector<bool> removed_triangle;

    //the barycenter
    vector<double> coordsOfB;

    void find_adjacencies();
    void find_incidencies();
    void link_adj (int t1, int t2);
    int valid_vertex(int v);
    int valid_triangle(int t);
};

template<class V, class T>
Mesh<V,T>::Mesh() {
    vertices = vector<V>();
    topSimplexes = vector<T>();
}

template<class V, class T>
Mesh<V,T>::Mesh(const Mesh& orig) {
    this->vertices = orig.vertices;
    this->topSimplexes = orig.topSimplexes;
}

template<class V, class T>
Mesh<V,T>::~Mesh() {
    vertices.clear();
    topSimplexes.clear();
}

template<class V, class T>
V& Mesh<V,T>::getVertex(int id){
    return this->vertices.at(id);
}

template<class V, class T>
int Mesh<V,T>::getNumVertex(){
    return this->vertices.size();
}

template<class V, class T>
void Mesh<V,T>::addVertex(V& v){
    this->vertices.push_back(v);
}

template<class V, class T>
T& Mesh<V,T>::getTopSimplex(int id){
    return this->topSimplexes.at(id);
}

template<class V, class T>
int Mesh<V,T>::getTopSimplexesNum(){
    return this->topSimplexes.size();
}

template<class V, class T>
void Mesh<V,T>::reserveVectorSpace(int numV, int numT){
    this->vertices.reserve(numV);
    this->topSimplexes.reserve(numT);
}

template<class V, class T>
void Mesh<V,T>::reserveVectorSpace_Vertices(int numV)
{
    this->vertices.reserve(numV);
}

template<class V, class T>
void Mesh<V,T>::reserveVectorSpace_TopSimplexes(int numT)
{
    this->topSimplexes.reserve(numT);
}

template<class V, class T>
void Mesh<V,T>::addTopSimplex(T& t){
    this->topSimplexes.push_back(t);
}

template<class V, class T>
void Mesh<V,T>::build()
{
    this->find_adjacencies();
    this->find_incidencies();
    removed_triangle = vector<bool>(getTopSimplexesNum(),false);
    removed_vertex = vector<bool>(getNumVertex(),false);
}

template<class V, class T>
void Mesh<V,T>::find_adjacencies()
{
    aux *tr_vec ;
    int i, j, k;
    int t1, t2;
    int v1, v2;

    tr_vec = (aux*) calloc( this->getTopSimplexesNum()*3 , sizeof(aux) ) ;
    if (!tr_vec)
    {
        cerr << "malloc fallita in find_adj" <<endl;
        return;
    }

    k = 0;

    for (j=0; j<this->getTopSimplexesNum(); j++)
    {
        for (i=0;i<3;i++)
        {
            tr_vec[k].t = j;
            this->getTopSimplex(j).setTT(i,-1);
            //            tl->elem[j]->adj[i] = -1;
            v1 = this->getTopSimplex(j).TV(i);
            //            v1 = tl->elem[j]->v[i];
            v2 = this->getTopSimplex(j).TV((i+1)%3);
            //            v2 = tl->elem[j]->v[(i+1)%3];
            if (v1<v2) {  tr_vec[k].v1 = v1; tr_vec[k].v2 = v2;  }
            else {  tr_vec[k].v1 = v2; tr_vec[k].v2 = v1;  }
            k++;
        }
    }

    qsort(tr_vec,3*this->getTopSimplexesNum(),sizeof(aux),cmp_aux) ;

    for(k=0;k<3*this->getTopSimplexesNum()-1;k++)
    {
        if ( is_endpoint(tr_vec[k].v1,tr_vec[k+1].v1,tr_vec[k+1].v2) &&
             is_endpoint(tr_vec[k].v2,tr_vec[k+1].v1,tr_vec[k+1].v2) )
        {  /* i due triangoli hanno lo stesso lato */
            t1 = tr_vec[k].t;
            t2 = tr_vec[k+1].t;
            link_adj(t1,t2);
        }
    }
    free(tr_vec) ;
    return;
}

template<class V, class T>
void Mesh<V,T>::find_incidencies()
{
    int i, t;

    for (t=0;t<this->getTopSimplexesNum();t++)
    {
        for (i=0;i<3;i++)
        {
            if(this->getVertex(this->getTopSimplex(t).TV(i)).VTstar()==-1)
                this->getVertex(this->getTopSimplex(t).TV(i)).VTstar(t);
            else if(this->getTopSimplex(t).TT((i+2)%3)==-1)
                this->getVertex(this->getTopSimplex(t).TV(i)).VTstar(t);
        }
    }
    return;
}

template<class V, class T>
void Mesh<V,T>::link_adj(int t1, int t2)
        /* Lega t1 come adiacente di t2 e viceversa */
{
    int i, j, k, pos1[2], pos2[2];
    if (valid_triangle(t1) && valid_triangle(t2))
    {
        k = 0;
        for (i=0; ((i<3)&&(k<2)); i++)
        {
            for (j=0; ((j<3)&&(k<2)); j++)
            {
                if(this->getTopSimplex(t1).TV(i) == this->getTopSimplex(t2).TV(j))
                {
                    pos1[k] = i;
                    pos2[k] = j;
                    k++;
                }
            }
        }
        if (k<2)
        {
            cerr << "error in link_adj" <<endl;
        }
        else
        {
            this->getTopSimplex(t1).setTT(3-pos1[0]-pos1[1],t2);
            this->getTopSimplex(t2).setTT(3-pos2[0]-pos2[1],t1);
        }
    }
    else
        cerr << "Error in link_adj: almeno uno dei triangoli e' nullo" << endl;
}

template<class V, class T>
int Mesh<V,T>::valid_vertex(int v)
{
    return ( (v>=0) && (v<this->getNumVertex()) );
}

template<class V, class T>
int Mesh<V,T>::valid_triangle(int t)
{
    return ( (t>=0) && (t<this->getTopSimplexesNum()) );
}

template<class V, class T>
vector<int> Mesh<V,T>::VT(int center)
{
    vector<int> triangles;
    int pred = -1;
    int current = this->getVertex(center).VTstar();
    triangles.push_back(this->getVertex(center).VTstar());

    int k=-1;

    //cerco la posizione del vertice nell'array dei vertici del triangolo
    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }

    //scelgo un giro a caso da prendere
    pred = current;
    current = this->getTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;
    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            isBorder = true;
            break;
        }
        else
            triangles.push_back(current);

        k=-1;
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }        

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    //se sono in un bordo ciclo anche nell'altro senso per recuperare i triangoli mancanti
    if(isBorder)
    {
        pred = this->getVertex(center).VTstar();
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(pred).getVerticesNum();i++)
        {
            if(this->getTopSimplex(pred).TV(i) == center)
            {
                k = i;
                break;
            }
        }
        current = this->getTopSimplex(pred).TT((k+2)%3);

        while(1)
        {
            if(current == -1)
                break;
            else
                triangles.push_back(current);

            k=-1;
            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
            {
                if(this->getTopSimplex(current).TV(i) == center)
                {
                    k = i;
                    break;
                }
            }


            if(this->getTopSimplex(current).TT((k+1)%3) == pred)
            {
                pred = current;
                current = this->getTopSimplex(current).TT((k+2)%3);
            }
            else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
            {
                pred = current;
                current = this->getTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return triangles;
}


template<class V, class T>
bool Mesh<V,T>::isBoundary(int center)
{
    vector<int> triangles;
    int pred = -1;
    int current = this->getVertex(center).VTstar();
    triangles.push_back(this->getVertex(center).VTstar());

    int k=-1;

    //cerco la posizione del vertice nell'array dei vertici del triangolo
    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }

    //scelgo un giro a caso da prendere
    pred = current;
    current = this->getTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;
    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            return true;
        }
        else
            triangles.push_back(current);

        k=-1;
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    //se sono in un bordo ciclo anche nell'altro senso per recuperare i triangoli mancanti
    return false;
}

template<class V, class T>
vector<Edge*> Mesh<V,T>::VE(int center)
{
    vector<Edge*> edges;
    int pred = -1;
    int current = this->getVertex(center).VTstar();

    int k=-1;
    //cerco la posizione del vertice nell'array dei vertici del triangolo
    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }
    edges.push_back(this->getTopSimplex(current).TE((k+1)%3));

    //scelgo un giro a caso da prendere
    pred = current;
    current = this->getTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;

    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            isBorder = true;
            break;
        }

        k=-1;
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            edges.push_back(this->getTopSimplex(current).TE((k+2)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            edges.push_back(this->getTopSimplex(current).TE((k+1)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    //se sono in un bordo ciclo anche nell'altro senso per recuperare i triangoli mancanti
    if(isBorder)
    {
        pred = this->getVertex(center).VTstar();
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(pred).getVerticesNum();i++)
        {
            if(this->getTopSimplex(pred).TV(i) == center)
            {
                k = i;
                break;
            }
        }
        edges.push_back(this->getTopSimplex(pred).TE((k+2)%3));
        current = this->getTopSimplex(pred).TT((k+2)%3);

        while(1)
        {
            if(current == -1)
                break;

            k=-1;
            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
            {
                if(this->getTopSimplex(current).TV(i) == center)
                {
                    k = i;
                    break;
                }
            }            

            if(this->getTopSimplex(current).TT((k+1)%3) == pred)
            {
                edges.push_back(this->getTopSimplex(current).TE((k+2)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+2)%3);
            }
            else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
            {
                edges.push_back(this->getTopSimplex(current).TE((k+1)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return edges;
}

template<class V, class T>
vector<int> Mesh<V,T>::VV(int center)
{
    vector<int> vertices;
    int pred = -1;
    int current = this->getVertex(center).VTstar();

    int k=-1;
    //cerco la posizione del vertice nell'array dei vertici del triangolo
    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
    {
        if(this->getTopSimplex(current).TV(i) == center)
        {
            k = i;
            break;
        }
    }    
    vertices.push_back(this->getTopSimplex(current).TV((k+1)%3));

    //scelgo un giro a caso da prendere
    pred = current;
    current = this->getTopSimplex(current).TT((k+2)%3);

    bool isBorder = false;

    while(1)
    {
        if(current == this->getVertex(center).VTstar())
            break;
        else if(current == -1)
        {
            isBorder = true;
            break;
        }

        k=-1;
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
        {
            if(this->getTopSimplex(current).TV(i) == center)
            {
                k = i;
                break;
            }
        }        

        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
        {
            vertices.push_back(this->getTopSimplex(current).TV((k+1)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+2)%3);
        }
        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
        {
            vertices.push_back(this->getTopSimplex(current).TV((k+2)%3));
            pred = current;
            current = this->getTopSimplex(current).TT((k+1)%3);
        }
    }

    //se sono in un bordo ciclo anche nell'altro senso per recuperare i triangoli mancanti
    if(isBorder)
    {
        pred = this->getVertex(center).VTstar();
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<this->getTopSimplex(pred).getVerticesNum();i++)
        {
            if(this->getTopSimplex(pred).TV(i) == center)
            {
                k = i;
                break;
            }
        }
        vertices.push_back(this->getTopSimplex(pred).TV((k+2)%3));
        current = this->getTopSimplex(pred).TT((k+1)%3);

        while(1)
        {
            if(current == -1)
                break;

            k=-1;
            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
            {
                if(this->getTopSimplex(current).TV(i) == center)
                {
                    k = i;
                    break;
                }
            }            

            if(this->getTopSimplex(current).TT((k+1)%3) == pred)
            {
                vertices.push_back(this->getTopSimplex(current).TV((k+1)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+2)%3);
            }
            else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
            {
                vertices.push_back(this->getTopSimplex(current).TV((k+2)%3));
                pred = current;
                current = this->getTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return vertices;
}

template<class V, class T>
vector<int> Mesh<V,T>::ET(Edge e)
{
    vector<int> triangles;
    vector<int> vtcomplete_triangles = this->VT(e.EV(0));
    for(unsigned int i=0;i<vtcomplete_triangles.size();i++)
    {
        for(int j=0;j<this->getTopSimplex(vtcomplete_triangles.at(i)).getVerticesNum();j++)
        {
            if(e.EV(1) == this->getTopSimplex(vtcomplete_triangles.at(i)).TV(j))
            {
                triangles.push_back(vtcomplete_triangles.at(i));
                break;
            }
        }
    }

    return triangles;
}

template<class V, class T>
vector<Edge*> Mesh<V,T>::EE(Edge e)
{
    vector<Edge*> edges;
    vector<Edge*> edge0 = this->VE(e.EV(0));
    vector<Edge*> edge1 = this->VE(e.EV(1));
    for(unsigned int i=0;i<edge0.size();i++)
    {
        if(e != *edge0.at(i))
            edges.push_back(edge0.at(i));
    }
    for(unsigned int i=0;i<edge1.size();i++)
    {
        if(e != *edge1.at(i))
            edges.push_back(edge1.at(i));
    }

    return edges;
}

template<class V, class T>
template<class C> bool Mesh<V,T> :: isIntoVector(C c, vector<C> &c_vect)
{
    for(unsigned int i=0; i<c_vect.size();i++)
        if(c == c_vect.at(i))
            return true;
    return false;
}


template<class V, class T>
void Mesh<V,T> :: reorder_triangulation(){

    vector<V> new_vertices;
    vector<T> new_triangle;

    vector<int> new_vertex_index = vector<int>(getNumVertex(),-1);
    vector<int> new_triangle_index = vector<int>(getTopSimplexesNum() ,-1);

    int verticesNum=0;
    int trianglesNum=0;
    for(int i=0; i<getNumVertex(); i++){

        if(!removed_vertex[i]){
            new_vertices.push_back(getVertex(i));
            assert(is_alive(getVertex(i).VTstar()));
            new_vertex_index[i]=verticesNum++;
        }
    }

    for(int i=0; i<getTopSimplexesNum(); i++){

        if(!removed_triangle[i]){
            for(int j=0; j<3; j++){
                assert(getTopSimplex(i).TT(j)==-1 || is_alive(getTopSimplex(i).TT(j)));
                assert(!removed_vertex[getTopSimplex(i).TV(j)]);

                getTopSimplex(i).setTV(j,new_vertex_index[getTopSimplex(i).TV(j)]);
            }
            new_triangle.push_back(getTopSimplex(i));
            new_triangle_index[i]=trianglesNum++;
        }
    }



    cout << "triangoli prima " << topSimplexes.size() << " dopo:" << new_triangle.size() << endl;
    cout << "vertices prima " << vertices.size() << " dopo:" << new_vertices.size() << endl;

    vertices = new_vertices;
    topSimplexes = new_triangle;

    removed_triangle = vector<bool>(new_triangle.size(), false);
    removed_vertex = vector<bool>(new_vertices.size(), false);

    for(int i=0; i<getNumVertex(); i++){
        getVertex(i).VTstar(new_triangle_index[getVertex(i).VTstar()]);
    }

    for(int i=0; i<getTopSimplexesNum(); i++){

            for(int j=0; j<3; j++){
                if(getTopSimplex(i).TT(j) != -1)
                    getTopSimplex(i).setTT(j, new_triangle_index[getTopSimplex(i).TT(j)]);
            }

    }

}


//return the index of vertex v in triangle t
template< class V, class T> int Mesh<V,T>::VIndexInT(int v, int t)
{
    for(int i=0; i<3; i++)
        if (this->getTopSimplex(t).TV(i)==v) return i;
    return -1;
}

//compute cosinus of angle formed by 3 vertices
template<class V, class T> double Mesh<V,T>::cosAngle(int v1, int v2, int v3)
{
    double norm21=this->getVertex(v2).norma(this->getVertex(v1));
    double norm23=this->getVertex(v2).norma(this->getVertex(v3));
    double product=this->getVertex(v2).prodscal(this->getVertex(v1),this->getVertex(v3));
    double costeta=product/(norm21*norm23);
    if (costeta>1.0) costeta = 1.0;
    if (costeta<-1.0) costeta = -1.0;
    return costeta;
}

//compute total area of triangles incident in v
template<class V, class T> double Mesh<V,T>::FanArea(int v)
{
    int t;
    double a = 0.0;

    t = this->getVertex(v).VTstar();
    if (t!=-1) do
    {
        a += TArea(t);
        t = NextTAroundV(t,v,CCW);
    }
    while ( (t!=-1) && (t!=this->getVertex(v).VTstar()) );


    if(t==-1)
    {
        t = this->getVertex(v).VTstar();
        if (t!=-1)
        {
            //non calcolo l'area per vtstar perchè è già stato fatto prima
            t = NextTAroundV(t,v,CW);
        }
        if (t!=-1) do
        {
            a += TArea(t);
            t = NextTAroundV(t,v,CW);
        }
        while ( (t!=-1) && (t!=this->getVertex(v).VTstar()) );
    }
    return a;
}

//compute triangle area
template<class V, class T> double Mesh<V,T>::TArea(int t)
{
    int v1,v2,v3;
    double prodscaluv;
    double normau;
    double normav;
    double cosalpha;
    double senalpha;

    v1=this->getTopSimplex(t).TV(0);
    v2=this->getTopSimplex(t).TV(1);
    v3=this->getTopSimplex(t).TV(2);
    prodscaluv=this->getVertex(v1).prodscal(this->getVertex(v2),this->getVertex(v3));
    normau=this->getVertex(v1).Norm(this->getVertex(v2));
    normav=this->getVertex(v1).Norm(this->getVertex(v3));
    cosalpha = prodscaluv / (normau*normav);
    senalpha=sqrt(1-(cosalpha*cosalpha));
    if(cosalpha*cosalpha > 1) senalpha = 0.0001;
    //if(isnan(senalpha)) senalpha=0.0001;

    return (normau*normav*senalpha)/2;
}

//returns the area of the whole mesh
template<class V, class T> double Mesh<V,T>::MArea(){

    double accumulator = 0.0;

    for(unsigned int k = 0; k < this->getTopSimplexesNum(); k++){
        if(this->TArea(k)>=0)
            accumulator += this->TArea(k);
    }
    return accumulator;
}

//returns the barycenter of the mesh
template<class V, class T> void Mesh<V,T>::Barycenter(){

    double xx=0.0, yy=0.0, zz=0.0;

    for(unsigned int kk = 0; kk < this->getNumVertex(); kk++){
        xx += this->getVertex(kk).getX();
        yy += this->getVertex(kk).getY();
        zz += this->getVertex(kk).getZ();
    }

    coordsOfB.push_back(xx/this->getNumVertex());
    coordsOfB.push_back(yy/this->getNumVertex());
    coordsOfB.push_back(zz/this->getNumVertex());
}

//return the next triangle incident in vertex v, starting form
//triangle t and moving in sense 'verso'
template<class V, class T> int Mesh<V,T>::NextTAroundV(int t, int v, versus verso)
{
    T& tr = this->getTopSimplex(t);
    int pos = VIndexInT(v, t);
    if(verso==CW)
        return tr.TT((pos+1)%3);
    else
        return tr.TT((pos+2)%3);
}

//compute barycentric area of triangle t centered at
//vertex v, where v is a vertex of t.
//compute Voronoi-barycentric area around vertex v
template<class V, class T> double Mesh<V,T>::VoronoiBarycentricArea(int v)
{
 double a = 0.0;
 int t = this->getVertex(v).VTstar();
 if (t!=-1) do
 {
   a += VoronoiBarycentricArea(v, t);
   t = NextTAroundV(t,v,CCW);
 }
 while ( (t!=-1) && (t!=this->getVertex(v).VTstar()) );
 return a;
}

//compute the contribution of triangle t to the mixed area around
//vertex v, where v is a vertex of t.
//the contribution of t is the Voronoi area if
//all angles are acute, or the barycentric area otherwise
template<class V, class T> double Mesh<V,T>::VoronoiBarycentricArea(int v, int t)
{
   int i, v1, v2;
   double cos_a, cos_b, cos_c;

// triangle t is v v1 v2
// a,b,c are the angles in v, v1, v2
   i = VIndexInT(v,t);
   v1 = this->getTopSimplex(t).TV((i+1)%3);
   v2 = this->getTopSimplex(t).TV((i+2)%3);
   cos_a = cosAngle(v2,v,v1);
   cos_b = cosAngle(v,v1,v2);
   cos_c = cosAngle(v1,v2,v);

   if ((cos_a<0.0)||(cos_b<0.0)||(cos_c<0.0)) // one is obtuse
   {
     if (cos_a<0.0) // the obtuse angle is in v
       return 0.5 * TArea(t);
     else return 0.25 * TArea(t);
   }
   else // Voronoi area
   {
     double sin_b, cot_b, sin_c, cot_c, e1, e2;
     /* if angles are zero area is zero */
     if (cos_b==1.0) return 0.0;
     if (cos_c==1.0) return 0.0;

     sin_b = sin(acos(cos_b));
     sin_c = sin(acos(cos_c));

     cot_b = cos_b / sin_b;
     cot_c = cos_c / sin_c;

     e1 = this->getVertex(v1).norma(this->getVertex(v));
     e2 = this->getVertex(v2).norma(this->getVertex(v));
     return (cot_c*e1*e1 + cot_b*e2*e2) / 8.0;
  }
}

//template<class V, class T>
//bool Mesh<V,T>::isBoundary(int center)
//{
//    vector<int> triangles;
//    int pred = -1;
//    int current = this->getVertex(center).VTstar();
//    triangles.push_back(this->getVertex(center).VTstar());

//    int k=-1;

//    //cerco la posizione del vertice nell'array dei vertici del triangolo
//    for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
//    {
//        if(this->getTopSimplex(current).TV(i) == center)
//        {
//            k = i;
//            break;
//        }
//    }

//    //scelgo un giro a caso da prendere
//    pred = current;
//    current = this->getTopSimplex(current).TT((k+1)%3);

////    bool isBorder = false;
//    while(1)
//    {
//        if(current == this->getVertex(center).VTstar())
//            break;
//        else if(current == -1)
//        {
//            return true;
//        }
//        else
//            triangles.push_back(current);

//        k=-1;
//        //cerco la posizione del vertice nell'array dei vertici del triangolo
//        for(int i=0;i<this->getTopSimplex(current).getVerticesNum();i++)
//        {
//            if(this->getTopSimplex(current).TV(i) == center)
//            {
//                k = i;
//                break;
//            }
//        }

//        if(this->getTopSimplex(current).TT((k+1)%3) == pred)
//        {
//            pred = current;
//            current = this->getTopSimplex(current).TT((k+2)%3);
//        }
//        else if(this->getTopSimplex(current).TT((k+2)%3) == pred)
//        {
//            pred = current;
//            current = this->getTopSimplex(current).TT((k+1)%3);
//        }
//    }

//    //se sono in un bordo ciclo anche nell'altro senso per recuperare i triangoli mancanti
//    return false;
//}


#endif // MESH_H
