#ifndef MESH_H
#define MESH_H

#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <assert.h>

#include "Sorting.h"
#include "Edge.h"
#include "Vertex3D.h"

/*
True iff v is an endpoint of edge v1-v2
*/
#define is_endpoint(v,v1,v2) ( ( (v)==(v1) ) || ( (v)==(v2) ) )

enum versus { CW=0 , CCW=1 };

/// A class representing a generic mesh parametrized by the type of top simplexes
template<class T>
class Mesh
{
public:
    ///A constructor method
    Mesh();
    ///A constructor method
    Mesh(const Mesh& pOrig);
    ///A destructor method
    virtual ~Mesh();
    ///A public method that returns the vertex at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a Vertex&, the vertex at the id-th position in the list
     */
    inline Vertex3D& GetVertex(int pId) { return m_Vertices[pId]; }
    ///A public method that returns the number of mesh vertices
    /*!
     * \return an integer, representing the number of vertices
     */
    inline int GetNumVertices() { return m_Vertices.size(); }

    ///A public method that adds a vertex to the vertices list
    /*!
     * \param v a Vertex& argument, representing the vertex to add
     */
    inline void AddVertex(const Vertex3D& pV) { m_Vertices.push_back(pV); }

    ///A public method that returns the face at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a T&, the face at the id-th position in the list
     */
    inline T& GetTopSimplex(int pId) { return m_TopSimplexes[pId]; }
    ///A public method that returns the number of top simplexes
    /*!
     * \return an integer, representing the number of top simplexes
     */
    inline int GetNumberOfTopSimplexes() { return m_TopSimplexes.size(); }
    inline int GetNumberOfTopSimplexes() const { return m_TopSimplexes.size(); }
    ///A public method that adds a top simplex to the top simplexes list
    /*!
     * \param t a T& argument, representing the top simplex to add
     */
    inline void AddTopSimplex(const T& t) { m_TopSimplexes.push_back(t); }

    ///A public method that initializes the space needed by the vertices and top simplexes lists
    /*!
     * \param numV an integer, represents the number of vertices
     * \param numT an integer, represents the number of top simplexes
     */
    void ReserveVectorSpace(int pNumV, int pNumT);
    void ReserveVectorSpace_TopSimplexes(int pNumT);
    void ReserveVectorSpace_Vertices(int pNumV);

    void Build();

    bool HalfedgeCollapse(int pV1, int pV2, int pT1, int pT2);
    bool ConvexNeighborhood(int pV1, int pV2, int pT1, int pT2);


    std::vector<int> VT(int pCenter);
    std::vector<Edge*> VE(int pCenter);
    std::vector<int> VV(int pCenter);

    bool IsBoundary(int pCenter);
    inline bool IsAlive(int pF){return !m_RemovedTriangle[pF];}
    inline bool IsVAlive(int pV){return !m_RemovedVertex[pV];}
    void ReorderTriangulation();


    //compute cosinus of angle formed by 3 vertices
    double CosAngle(int pV1, int pV2, int pV3);
    //compute total area of triangles incident in v
    double FanArea(int v);
    //compute triangle area
    double TriangleArea(int t);
    //computes area of the whole mesh
    double MeshArea();
    //barycenter
    void Barycenter();

    inline glm::vec3 GetBarycenter() const { return m_Barycenter; }
    float GetGoundingBoxDiagonal();
    Vertex3D HalfPoint(unsigned int pIndexV1, unsigned int pIndexV2);

    int VIndexInT(int pV, int pT);
    int NextTAroundV(int pT, int pV, versus pVerso);

    //compute the mixed (Voronoi-barycentric) area around vertex v
    double VoronoiBarycentricArea(int pV);
    //compute either the Voronoi area (if all angles of t are acute) or the
    //barycentric area (otherwise) of triangle t centered in its vertex v
    double VoronoiBarycentricArea(int pV, int pT);


    std::vector<int> ET(const Edge& pE);
    std::vector<Edge*> EE(const Edge& pE);


    //Utility Function
    template<class C> bool IsIntoVector(C c, const std::vector<C>& c_vect);


private:
    void FindAdjacencies();
    void FindIncidencies();
    void LinkAdjacencies(int pT1, int pT2);
    int ValidVertex(int pV);
    int ValidTriangle(int pT);

private:
    ///A private varible representing the vertices list of the mesh
    std::vector<Vertex3D> m_Vertices;
    ///A private varible representing the top simplexes list of the mesh
    std::vector<T> m_TopSimplexes;

    std::vector<bool> m_RemovedVertex;
    std::vector<bool> m_RemovedTriangle;

    //the barycenter
    glm::vec3 m_Barycenter;

};


template<class T>
Mesh<T>::Mesh()
    :
    m_Vertices(std::vector<Vertex3D>())
    , m_TopSimplexes(std::vector<T>())
{

}

template<class T>
Mesh<T>::Mesh(const Mesh& pOrig)
    :
    m_Vertices(pOrig.m_Vertices)
    , m_TopSimplexes(pOrig.m_TopSimplexes)
{
}

template<class T>
Mesh<T>::~Mesh()
{
    m_Vertices.clear();
    m_TopSimplexes.clear();
}

template<class T>
void Mesh<T>::ReserveVectorSpace(int numV, int numT)
{
    m_Vertices.reserve(numV);
    m_TopSimplexes.reserve(numT);
}

template<class T>
void Mesh<T>::ReserveVectorSpace_Vertices(int numV)
{
    m_Vertices.reserve(numV);
}

template<class T>
void Mesh<T>::ReserveVectorSpace_TopSimplexes(int numT)
{
    m_TopSimplexes.reserve(numT);
}


template< class T>
void Mesh<T>::Build()
{
    FindAdjacencies();
    FindIncidencies();
    m_RemovedTriangle = std::vector<bool>(GetNumberOfTopSimplexes(),false);
    m_RemovedVertex = std::vector<bool>(GetNumVertices(),false);
}

template<class T>
void Mesh<T>::FindAdjacencies()
{
    //aux *tr_vec ;
    std::vector<aux> tr_vec;

    tr_vec.reserve(GetNumberOfTopSimplexes()*3);
    int k = 0;

    for (int j=0; j<GetNumberOfTopSimplexes(); j++)
    {
        for (int i=0;i<3;i++)
        {
            tr_vec[k].t = j;
            GetTopSimplex(j).setTT(i,-1);
            //            tl->elem[j]->adj[i] = -1;
            int v1 = GetTopSimplex(j).TV(i);
            //            v1 = tl->elem[j]->v[i];
            int v2 = GetTopSimplex(j).TV((i+1)%3);
            //            v2 = tl->elem[j]->v[(i+1)%3];
            if (v1<v2)
            {
                tr_vec[k].v1 = v1;
                tr_vec[k].v2 = v2;
            }
            else
            {
                tr_vec[k].v1 = v2;
                tr_vec[k].v2 = v1;
            }
            k++;
        }
    }

    std::sort(tr_vec.begin(), tr_vec.end(), CompareAux);

    for(k=0; k<3*GetNumberOfTopSimplexes()-1; k++)
    {
        if ( is_endpoint(tr_vec[k].v1,tr_vec[k+1].v1,tr_vec[k+1].v2) &&
             is_endpoint(tr_vec[k].v2,tr_vec[k+1].v1,tr_vec[k+1].v2) )
        {  /* i due triangoli hanno lo stesso lato */
            int t1 = tr_vec[k].t;
            int t2 = tr_vec[k+1].t;
            LinkAdjacencies(t1,t2);
        }
    }
    tr_vec.clear();
    return;
}

template<class T>
void Mesh<T>::FindIncidencies()
{
    for (int t=0; t<GetNumberOfTopSimplexes(); t++)
    {
        for (int i=0; i<3; i++)
        {
            if(GetVertex(GetTopSimplex(t).TV(i)).VTstar()==-1)
                GetVertex(GetTopSimplex(t).TV(i)).VTstar(t);
            else if(GetTopSimplex(t).TT((i+2)%3)==-1)
                GetVertex(GetTopSimplex(t).TV(i)).VTstar(t);
        }
    }
    return;
}

template<class T>
void Mesh<T>::LinkAdjacencies(int pT1, int pT2)
        /* Lega t1 come adiacente di t2 e viceversa */
{
    if (ValidTriangle(pT1) && ValidTriangle(pT2))
    {
        int k = 0;
        int pos1[2], pos2[2];
        for (int i=0; ((i<3)&&(k<2)); i++)
        {
            for (int j=0; ((j<3)&&(k<2)); j++)
            {
                if(GetTopSimplex(pT1).TV(i) == GetTopSimplex(pT2).TV(j))
                {
                    pos1[k] = i;
                    pos2[k] = j;
                    k++;
                }
            }
        }

        if (k<2)
        {
            std::cerr << "error in link_adj" << std::endl;
        }
        else
        {
            GetTopSimplex(pT1).setTT(3-pos1[0]-pos1[1],pT2);
            GetTopSimplex(pT2).setTT(3-pos2[0]-pos2[1],pT1);
        }
    }
    else
        std::cerr << "Error in link_adj: at least one of the triangles is NULL" << std::endl;
}

template<class T>
int Mesh<T>::ValidVertex(int v)
{
    return ( (v>=0) && (v < GetNumVertices()) );
}

template<class T>
int Mesh<T>::ValidTriangle(int t)
{
    return ( (t>=0) && (t < GetNumberOfTopSimplexes()) );
}

template<class T>
std::vector<int> Mesh<T>::VT(int pCenter)
{
    std::vector<int> triangles;
    int pred = -1;
    int current = GetVertex(pCenter).VTstar();
    triangles.push_back(GetVertex(pCenter).VTstar());

    int k=-1;

    //look for the vertex position in the triangle vertices array
    for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
    {
        if(GetTopSimplex(current).TV(i) == pCenter)
        {
            k = i;
            break;
        }
    }

    //choose a random turn
    pred = current;
    current = GetTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;
    while(1)
    {
        if(current == GetVertex(pCenter).VTstar())
        {
            break;
        }
        else if(current == -1)
        {
            isBorder = true;
            break;
        }
        else
        {
            triangles.push_back(current);
        }

        k=-1;
        //look for the vertex position within the triangle vertices array
        for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
        {
            if(GetTopSimplex(current).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }

        if(GetTopSimplex(current).TT((k+1)%3) == pred)
        {
            pred = current;
            current = GetTopSimplex(current).TT((k+2)%3);
        }
        else if(GetTopSimplex(current).TT((k+2)%3) == pred)
        {
            pred = current;
            current = GetTopSimplex(current).TT((k+1)%3);
        }
    }

    //if I'm on the border I loop in the other direction to account for triangles I might have missed
    if(isBorder)
    {
        pred = GetVertex(pCenter).VTstar();
        //search the vertex position in the triangle vertices array
        for(int i=0;i<GetTopSimplex(pred).getVerticesNum();i++)
        {
            if(GetTopSimplex(pred).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }
        current = GetTopSimplex(pred).TT((k+2)%3);

        while(1)
        {
            if(current == -1)
                break;
            else
                triangles.push_back(current);

            k=-1;
            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
            {
                if(GetTopSimplex(current).TV(i) == pCenter)
                {
                    k = i;
                    break;
                }
            }


            if(GetTopSimplex(current).TT((k+1)%3) == pred)
            {
                pred = current;
                current = GetTopSimplex(current).TT((k+2)%3);
            }
            else if(GetTopSimplex(current).TT((k+2)%3) == pred)
            {
                pred = current;
                current = GetTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return triangles;
}


template<class T>
bool Mesh<T>::IsBoundary(int pCenter)
{
    std::vector<int> triangles;
    int pred = -1;
    int current = GetVertex(pCenter).VTstar();
    triangles.push_back(GetVertex(pCenter).VTstar());

    int k=-1;

    //look for the vertex position in the triangle vertices array
    for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
    {
        if(GetTopSimplex(current).TV(i) == pCenter)
        {
            k = i;
            break;
        }
    }

    //choose a random direction to follow
    pred = current;
    current = GetTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;
    while(1)
    {
        if(current == GetVertex(pCenter).VTstar())
        {
            break;
        }
        else if(current == -1)
        {
            return true;
        }
        else
        {
            triangles.push_back(current);
        }

        k=-1;

        //look for the vertex position in the triangle vertices array
        for(int i=0; i<GetTopSimplex(current).getVerticesNum(); i++)
        {
            if(GetTopSimplex(current).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }

        if(GetTopSimplex(current).TT((k+1)%3) == pred)
        {
            pred = current;
            current = GetTopSimplex(current).TT((k+2)%3);
        }
        else if(GetTopSimplex(current).TT((k+2)%3) == pred)
        {
            pred = current;
            current = GetTopSimplex(current).TT((k+1)%3);
        }
    }

    return false;
}

template<class T>
std::vector<Edge*> Mesh<T>::VE(int pCenter)
{
    std::vector<Edge*> edges;
    int pred = -1;
    int current = GetVertex(pCenter).VTstar();

    int k=-1;
    //look for the vertex position in the triangle vertices array
    for(int i=0; i<GetTopSimplex(current).getVerticesNum(); i++)
    {
        if(GetTopSimplex(current).TV(i) == pCenter)
        {
            k = i;
            break;
        }
    }
    edges.push_back(GetTopSimplex(current).TE((k+1)%3));

    //choose a direction to loop in
    pred = current;
    current = GetTopSimplex(current).TT((k+1)%3);

    bool isBorder = false;

    while(1)
    {
        if(current == GetVertex(pCenter).VTstar())
        {
            break;
        }
        else if(current == -1)
        {
            isBorder = true;
            break;
        }

        k=-1;

        //look for the vertex position in the triangle vertices array
        for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
        {
            if(GetTopSimplex(current).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }

        if(GetTopSimplex(current).TT((k+1)%3) == pred)
        {
            edges.push_back(GetTopSimplex(current).TE((k+2)%3));
            pred = current;
            current = GetTopSimplex(current).TT((k+2)%3);
        }
        else if(GetTopSimplex(current).TT((k+2)%3) == pred)
        {
            edges.push_back(GetTopSimplex(current).TE((k+1)%3));
            pred = current;
            current = GetTopSimplex(current).TT((k+1)%3);
        }
    }

    //if I'm on a boorder, I loop in the other direction to look for missing triangles
    if(isBorder)
    {
        pred = GetVertex(pCenter).VTstar();
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0;i<GetTopSimplex(pred).getVerticesNum();i++)
        {
            if(GetTopSimplex(pred).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }
        edges.push_back(GetTopSimplex(pred).TE((k+2)%3));
        current = GetTopSimplex(pred).TT((k+2)%3);

        while(1)
        {
            if(current == -1)
                break;

            k=-1;
            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
            {
                if(GetTopSimplex(current).TV(i) == pCenter)
                {
                    k = i;
                    break;
                }
            }

            if(GetTopSimplex(current).TT((k+1)%3) == pred)
            {
                edges.push_back(GetTopSimplex(current).TE((k+2)%3));
                pred = current;
                current = GetTopSimplex(current).TT((k+2)%3);
            }
            else if(GetTopSimplex(current).TT((k+2)%3) == pred)
            {
                edges.push_back(GetTopSimplex(current).TE((k+1)%3));
                pred = current;
                current = GetTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return edges;
}

template<class T>
std::vector<int> Mesh<T>::VV(int pCenter)
{
    std::vector<int> vertices;
    int pred = -1;
    int current = GetVertex(pCenter).VTstar();

    int k=-1;
    //look for the vertex position in the triangle vertices array
    for(int i=0; i<GetTopSimplex(current).getVerticesNum(); i++)
    {
        if(GetTopSimplex(current).TV(i) == pCenter)
        {
            k = i;
            break;
        }
    }
    vertices.push_back(GetTopSimplex(current).TV((k+1)%3));

    //loop in a random direction
    pred = current;
    current = GetTopSimplex(current).TT((k+2)%3);

    bool isBorder = false;

    while(1)
    {
        if(current == GetVertex(pCenter).VTstar())
        {
            break;
        }
        else if(current == -1)
        {
            isBorder = true;
            break;
        }

        k=-1;

        //look for the vertex position in the triangle vertices array
        for(int i=0;i<GetTopSimplex(current).getVerticesNum();i++)
        {
            if(GetTopSimplex(current).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }

        if(GetTopSimplex(current).TT((k+1)%3) == pred)
        {
            vertices.push_back(GetTopSimplex(current).TV((k+1)%3));
            pred = current;
            current = GetTopSimplex(current).TT((k+2)%3);
        }
        else if(GetTopSimplex(current).TT((k+2)%3) == pred)
        {
            vertices.push_back(GetTopSimplex(current).TV((k+2)%3));
            pred = current;
            current = GetTopSimplex(current).TT((k+1)%3);
        }
    }

    //if I'm on a border I loop in the other direction too to look for missing triangles
    if(isBorder)
    {
        pred = GetVertex(pCenter).VTstar();
        //cerco la posizione del vertice nell'array dei vertici del triangolo
        for(int i=0; i<GetTopSimplex(pred).getVerticesNum(); i++)
        {
            if(GetTopSimplex(pred).TV(i) == pCenter)
            {
                k = i;
                break;
            }
        }
        vertices.push_back(GetTopSimplex(pred).TV((k+2)%3));
        current = GetTopSimplex(pred).TT((k+1)%3);

        while(1)
        {
            if(current == -1)
            {
                break;
            }

            k=-1;
            //look for the vertex position in the triangle vertices array
            for(int i=0; i<GetTopSimplex(current).getVerticesNum(); i++)
            {
                if(GetTopSimplex(current).TV(i) == pCenter)
                {
                    k = i;
                    break;
                }
            }

            if(GetTopSimplex(current).TT((k+1)%3) == pred)
            {
                vertices.push_back(GetTopSimplex(current).TV((k+1)%3));
                pred = current;
                current = GetTopSimplex(current).TT((k+2)%3);
            }
            else if(GetTopSimplex(current).TT((k+2)%3) == pred)
            {
                vertices.push_back(GetTopSimplex(current).TV((k+2)%3));
                pred = current;
                current = GetTopSimplex(current).TT((k+1)%3);
            }
        }
    }

    return vertices;
}

template<class T>
std::vector<int> Mesh<T>::ET(const Edge& pE)
{
    std::vector<int> triangles;
    std::vector<int> vtcomplete_triangles = this->VT(pE.EV(0));
    for(unsigned int i=0; i<vtcomplete_triangles.size(); i++)
    {
        for(int j=0;j<GetTopSimplex(vtcomplete_triangles.at(i)).getVerticesNum();j++)
        {
            if(pE.EV(1) == GetTopSimplex(vtcomplete_triangles.at(i)).TV(j))
            {
                triangles.push_back(vtcomplete_triangles.at(i));
                break;
            }
        }
    }

    return triangles;
}

template<class T>
std::vector<Edge*> Mesh<T>::EE(const Edge& pE)
{
    std::vector<Edge*> edges;
    std::vector<Edge*> edge0 = this->VE(pE.EV(0));
    std::vector<Edge*> edge1 = this->VE(pE.EV(1));
    for(unsigned int i=0;i<edge0.size();i++)
    {
        if(pE != *edge0.at(i))
            edges.push_back(edge0.at(i));
    }
    for(unsigned int i=0;i<edge1.size();i++)
    {
        if(pE != *edge1.at(i))
            edges.push_back(edge1.at(i));
    }

    return edges;
}

template<class T>
template<class C> bool Mesh<T> :: IsIntoVector(C c, const std::vector<C> &c_vect)
{
    for(unsigned int i=0; i<c_vect.size();i++)
        if(c == c_vect.at(i))
            return true;
    return false;
}


template<class T>
void Mesh<T> :: ReorderTriangulation()
{

    std::vector<Vertex3D> new_vertices;
    std::vector<T> new_triangle;

    std::vector<int> new_vertex_index = vector<int>(GetNumVertices(),-1);
    std::vector<int> new_triangle_index = vector<int>(GetNumberOfTopSimplexes() ,-1);

    int verticesNum=0;
    int trianglesNum=0;
    for(int i=0; i<GetNumVertices(); i++)
    {

        if(!m_RemovedVertex[i])
        {
            new_vertices.push_back(GetVertex(i));
            assert(is_alive(GetVertex(i).VTstar()));
            new_vertex_index[i]=verticesNum++;
        }
    }

    for(int i=0; i<GetNumberOfTopSimplexes(); i++)
    {

        if(!m_RemovedTriangle[i])
        {
            for(int j=0; j<3; j++)
            {
                assert(GetTopSimplex(i).TT(j)==-1 || is_alive(GetTopSimplex(i).TT(j)));
                assert(!m_RemovedVertex[GetTopSimplex(i).TV(j)]);

                GetTopSimplex(i).setTV(j,new_vertex_index[GetTopSimplex(i).TV(j)]);
            }
            new_triangle.push_back(GetTopSimplex(i));
            new_triangle_index[i]=trianglesNum++;
        }
    }

    m_Vertices = new_vertices;
    m_TopSimplexes = new_triangle;

    m_RemovedTriangle = std::vector<bool>(new_triangle.size(), false);
    m_RemovedVertex = std::vector<bool>(new_vertices.size(), false);

    for(int i=0; i<GetNumVertices(); i++)
    {
        GetVertex(i).VTstar(new_triangle_index[GetVertex(i).VTstar()]);
    }

    for(int i=0; i<GetNumberOfTopSimplexes(); i++)
    {

            for(int j=0; j<3; j++){
                if(GetTopSimplex(i).TT(j) != -1)
                    GetTopSimplex(i).setTT(j, new_triangle_index[GetTopSimplex(i).TT(j)]);
            }

    }

}


//return the index of vertex v in triangle t
template<class T> int Mesh<T>::VIndexInT(int v, int t)
{
    for(int i=0; i<3; i++)
    {
        if (GetTopSimplex(t).TV(i)==v)
        {
            return i;
        }
    }
    return -1;
}

//compute cosinus of angle formed by 3 vertices. The angle is pV1pV2pV3
template<class T> double Mesh<T>::CosAngle(int pV1, int pV2, int pV3)
{
    //sides adjacent to the angle
    double norm21 = GetVertex(pV2).SquaredDistance(GetVertex(pV1));
    if(norm21 == 0.0f)
    {
        norm21 = 0.01;
    }
    double norm23 = GetVertex(pV2).SquaredDistance(GetVertex(pV3));
    if(norm23 == 0.0f)
    {
        norm23 = 0.01;
    }

    //opposite side wrt the angle
    double norm13 = GetVertex(pV3).SquaredDistance(GetVertex(pV1));
    if(norm13 == 0.0f)
    {
        norm13 = 0.01;
    }

    double numerator = norm21 + norm23 - norm13;
    double denominator = 2.0 * sqrt(norm21 * norm23); //sqrt(a*b) = sqrt(a)*sqrt(b)

    return numerator / denominator;
}

// //compute total area of triangles incident in v
template<class T> double Mesh<T>::FanArea(int v)
{
    double a = 0.0;

    int t = GetVertex(v).VTstar();

    if (t!=-1) do
    {
            a += TriangleArea(t);
        t = NextTAroundV(t,v,CCW);
    }
    while ( (t!=-1) && (t!=GetVertex(v).VTstar()) );


    if(t==-1)
    {
        t = GetVertex(v).VTstar();
        if (t!=-1)
        {
            //non calcolo l'area per vtstar perchè è già stato fatto prima
            t = NextTAroundV(t,v,CW);
        }
        if (t!=-1) do
        {
                a += TriangleArea(t);
            t = NextTAroundV(t,v,CW);
        }
        while ( (t!=-1) && (t!=GetVertex(v).VTstar()) );
    }
    return a;
}

//compute triangle area
template<class T> double Mesh<T>::TriangleArea(int pT)
{
    int v1=GetTopSimplex(pT).TV(0);
    int v2=GetTopSimplex(pT).TV(1);
    int v3=GetTopSimplex(pT).TV(2);

    double v1v2 = GetVertex(v1).distance(GetVertex(v2));
    double v2v3 = GetVertex(v2).distance(GetVertex(v3));
    double v3v1 = GetVertex(v3).distance(GetVertex(v1));

    double halfPerimeter = (v1v2 + v2v3 + v3v1) / 2.0;

    return sqrt(halfPerimeter * (halfPerimeter-v1v2) * (halfPerimeter-v2v3) * (halfPerimeter-v3v1));
}

//returns the area of the whole mesh
template<class T> double Mesh<T>::MeshArea()
{

    double accumulator = 0.0;

    for(unsigned int k = 0; k < GetNumberOfTopSimplexes(); k++)
    {
        if(this->TriangleArea(k)>=0)
        {
            accumulator += TriangleArea(k);
        }
    }
    return accumulator;
}

//returns the barycenter of the mesh
template<class T> void Mesh<T>::Barycenter()
{

    double xx=0.0, yy=0.0, zz=0.0;

    for(unsigned int kk = 0; kk < GetNumVertices(); kk++)
    {
        xx += this->getVertex(kk).getX();
        yy += this->getVertex(kk).getY();
        zz += this->getVertex(kk).getZ();
    }

    m_Barycenter = glm::vec3(xx/GetNumVertices(), yy/GetNumVertices(), zz/GetNumVertices());
}

//return the next triangle incident in vertex v, starting form
//triangle t and moving in sense 'verso'
template<class T> int Mesh<T>::NextTAroundV(int t, int v, versus verso)
{
    T& tr = GetTopSimplex(t);
    int pos = VIndexInT(v, t);

    return (verso == CW) ? tr.TT((pos+1)%3) : tr.TT((pos+2)%3);
}

//compute barycentric area of triangle t centered at
//vertex v, where v is a vertex of t.
//compute Voronoi-barycentric area around vertex v
template<class T> double Mesh<T>::VoronoiBarycentricArea(int v)
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
template<class T> double Mesh<T>::VoronoiBarycentricArea(int v, int t)
{
   // triangle t is v v1 v2
   // a,b,c are the angles in v, v1, v2
   int i = VIndexInT(v,t);
   int v1 = GetTopSimplex(t).TV((i+1)%3);
   int v2 = GetTopSimplex(t).TV((i+2)%3);
   double cos_a = CosAngle(v2,v,v1);
   double cos_b = CosAngle(v,v1,v2);
   double cos_c = CosAngle(v1,v2,v);

   if ((cos_a<0.0)||(cos_b<0.0)||(cos_c<0.0)) // one is obtuse
   {
     if (cos_a<0.0) // the obtuse angle is in v
       return 0.5 * TriangleArea(t);
     else
         return 0.25 * TriangleArea(t);
   }
   else // Voronoi area
   {
     /* if angles are zero area is zero */
     if (cos_b==1.0) return 0.0;
     if (cos_c==1.0) return 0.0;

     double sin_b = sin(acos(cos_b));
     double sin_c = sin(acos(cos_c));

     double cot_b = cos_b / sin_b;
     double cot_c = cos_c / sin_c;

     double e1 = GetVertex(v1).norma(GetVertex(v));
     double e2 = GetVertex(v2).norma(GetVertex(v));
     return (cot_c*e1*e1 + cot_b*e2*e2) / 8.0;
  }
}


template<class T>
Vertex3D Mesh<T>::HalfPoint(unsigned int pIndexV1, unsigned int pIndexV2)
{
    glm::vec3 coords1 = GetVertex(pIndexV1).GetCoordinates();
    glm::vec3 coords2 = GetVertex(pIndexV2).GetCoordinates();

    glm::vec3 hp = (coords1 + coords2) / 2.0f;

    return Vertex3D(hp);
}


template<class T>
float Mesh<T>::GetGoundingBoxDiagonal()
{
    //Minimum
    glm::vec3 minCoordinate = GetVertex(0).GetCoordinates();
    //Maximum
    glm::vec3 maxCoordinate = GetVertex(0).GetCoordinates();

    for(int ii=0; ii< GetNumVertices(); ii++)
    {

        Vertex3D vv = GetVertex(ii);

        //Update if necessary
        if(vv.getX() < minCoordinate.x)
            minCoordinate.x = vv.getX();
        if(vv.getY() < minCoordinate.y)
            minCoordinate.y = vv.getY();
        if(vv.getZ() < minCoordinate.z)
            minCoordinate.z = vv.getZ();

        if(vv.getX() > maxCoordinate.x)
            maxCoordinate.x = vv.getX();
        if(vv.getY() > maxCoordinate.y)
            maxCoordinate.y = vv.getY();
        if(vv.getZ() > maxCoordinate.z)
            maxCoordinate.z = vv.getZ();
    }

    //diagonal length
    return glm::distance(minCoordinate, maxCoordinate);
}


#endif // MESH_H
