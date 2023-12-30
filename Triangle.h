#ifndef Triangle_H
#define Triangle_H

#include "Edge.h"


/**
 * @brief The Triangle class stores information about the triangles which compose a mesh
 */
class Triangle
{
public:
    Triangle();
    Triangle(int v1, int v2, int v3);

    int TV(int pos);
    int getVerticesNum();

    Edge* TE(int pos);

    int TT(int pos);
    void setTT(int pos, int adjId);

    inline void setTV(int pos, int v){
        m_Vertices[pos]=v;
    }

    bool operator== (const Triangle &p) const;


    inline bool operator!= (const Triangle p)
    {
       return !((*this)==(p));
    }

    inline int vertex_index(int v){
        for(int i=0; i<3; i++){
            if(m_Vertices[i] == v) return i;
        }
        return -1;
    }

    inline bool contains(int v){
        for(int i=0; i<3; i++){
            if(m_Vertices[i] == v) return true;
        }
        return false;
    }

    inline glm::ivec3 Vertices() const
    {
        return m_Vertices;
    }

    inline glm::ivec3 Adjacencies() const
    {
        return m_Adjacencies;
    }


private:
    glm::ivec3 m_Vertices;
    ///
    glm::ivec3 m_Adjacencies;
};



#endif // Triangle_H
