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
{
    this->vertices[0] = v1;
    this->vertices[1] = v2;
    this->vertices[2] = v3;
}

/**
 * @brief Triangle::TV
 * @param pos index of the vertex
 * @return the incident vertex in the pos-th position
 */
int Triangle::TV(int pos)
{
    return this->vertices[pos];
}

/**
 * @brief Triangle::TE
 * @param pos index of the edge
 * @return incident edge in pos-th position
 */
Edge* Triangle::TE(int pos)
{
    return new Edge(vertices[(pos+1)%3],vertices[(pos+2)%3]);
}

/**
 * @brief Triangle::TT
 * @param pos index of the triangle
 * @return adjacent triangle in pos-th position
 */
int Triangle::TT(int pos)
{
    return this->adj[pos];
}

/**
 * @brief Triangle::setTT
 * @param pos position (0,1,or 2)
 * @param adjId (index of ajdacent triangle)
 */
void Triangle::setTT(int pos, int adjId)
{
    this->adj[pos]=adjId;
}

/**
 * @brief Triangle::getVerticesNum
 * @return number of vertices in a triangle (is supposed to be 3)
 */
int Triangle::getVerticesNum()
{
    return 3;
}
