#include "Edge.h"

Edge::Edge() {}

/**
 * @brief Edge::Edge Constructor
 * @param v1 index of the first vertex
 * @param v2 index of the second vertex
 */
Edge::Edge(int v1, int v2)
{
    this->vertices[0] = v1;
    this->vertices[1] = v2;
}

/**
 * @brief Edge::EV implementation of the edge-vertex relation
 * @param pos (values 0 or 1) which one of the two vertices we want
 * @return index of the corresponding vertex
 */
int Edge::EV(int pos)
{
    return this->vertices[pos];
}
