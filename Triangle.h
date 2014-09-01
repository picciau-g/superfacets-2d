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
        vertices[pos]=v;
    }

    inline bool operator== (const Triangle &p)
    {
        bool b[3];
        b[0] = false; b[1] = false; b[2] = false;
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(!b[j] && p.vertices[i]==this->vertices[j])
                {
                    b[j] = true;
                    break;
                }
            }
        }

        return b[0] && b[1] && b[2];
    }

    inline bool operator!= (const Triangle p)
    {
       return !((*this)==(p));
    }

    inline int vertex_index(int v){
        for(int i=0; i<3; i++){
            if(vertices[i] == v) return i;
        }
        return -1;
    }

    inline bool contains(int v){
        for(int i=0; i<3; i++){
            if(vertices[i] == v) return true;
        }
        return false;
    }

private:
    int vertices[3];
    ///
    int adj[3];
};



#endif // Triangle_H
