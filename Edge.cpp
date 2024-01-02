#include "Edge.h"

Edge::Edge() {}

/**
 * @brief Edge::Edge Constructor
 * @param v1 index of the first vertex
 * @param v2 index of the second vertex
 */
Edge::Edge(int v1, int v2)
    :
    m_vertIndices(glm::vec2(v1,v2))
{
}

Edge::Edge(const Edge& pOth)
    :
    m_vertIndices(pOth.m_vertIndices)
{

}

/**
 * @brief Edge::EV implementation of the edge-vertex relation
 * @param pos (values 0 or 1) which one of the two vertices we want
 * @return index of the corresponding vertex
 */
int Edge::EV(int pos)
{
    return (pos == 0) ? m_vertIndices.x : m_vertIndices.y;
}

/**
 * @brief Edge::EV implementation of the edge-vertex relation
 * @param pos (values 0 or 1) which one of the two vertices we want
 * @return index of the corresponding vertex
 */
int Edge::EV(int pos) const
{
    return (pos == 0) ? m_vertIndices.x : m_vertIndices.y;
}


bool Edge::operator==(const Edge &p) const
{
    return m_vertIndices == p.VertexIndices() ||
           (m_vertIndices.x == p.VertexIndices().y && m_vertIndices.y == p.VertexIndices().x);
}
