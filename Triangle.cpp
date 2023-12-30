#include "Triangle.h"


Triangle::Triangle()
{
}

/**
 * @brief Triangle::Triangle class constructor
 * @param v1 first vertex
 * @param v2 second vertex
 * @param v3 third vertex
 */
Triangle::Triangle(int v1, int v2, int v3)
    :
    m_Vertices(glm::ivec3(v1, v2, v3))
{
}

/**
 * @brief Triangle::TV
 * @param pos index of the vertex
 * @return the incident vertex in the pos-th position
 */
int Triangle::TV(int pos)
{
    return m_Vertices[pos];
}

/**
 * @brief Triangle::TE
 * @param pos index of the edge
 * @return incident edge in pos-th position
 */
Edge* Triangle::TE(int pos)
{
    return new Edge(m_Vertices[(pos+1)%3],m_Vertices[(pos+2)%3]);
}

/**
 * @brief Triangle::TT
 * @param pos index of the triangle
 * @return adjacent triangle in pos-th position
 */
int Triangle::TT(int pos)
{
    return this->m_Adjacencies[pos];
}

/**
 * @brief Triangle::setTT
 * @param pos position (0,1,or 2)
 * @param adjId (index of ajdacent triangle)
 */
void Triangle::setTT(int pos, int adjId)
{
    this->m_Adjacencies[pos]=adjId;
}

/**
 * @brief Triangle::getVerticesNum
 * @return number of vertices in a triangle (is supposed to be 3)
 */
int Triangle::getVerticesNum()
{
    return 3;
}


bool Triangle::operator==(const Triangle& pOth) const
{
    return m_Vertices == pOth.Vertices() ||
           (m_Vertices.x == pOth.Vertices().y && m_Vertices.y == pOth.Vertices().z && m_Vertices.z == pOth.Vertices().x) ||
           (m_Vertices.x == pOth.Vertices().z && m_Vertices.z == pOth.Vertices().y && m_Vertices.y == pOth.Vertices().x);
}
