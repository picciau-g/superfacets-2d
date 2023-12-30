#ifndef EDGE_H
#define EDGE_H

/**
 * @brief The Edge class stores minimal information about edges of a mesh
 */

#include <glm/glm.hpp>

class Edge
{
public:
    Edge();
    Edge(int v1, int v2);
    Edge(const Edge& pOth);

    int EV(int pos);

    bool operator== (const Edge &p) const;

    inline glm::vec2 VertexIndices() const
    {
        return m_vertIndices;
    }

    inline bool operator!= (Edge &p) const
    {
       return !(p==(*this));
    }

private:
    glm::vec2 m_vertIndices;
};


#endif // EDGE_H
