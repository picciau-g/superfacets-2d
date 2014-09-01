#include "Edge.h"

Edge::Edge() {}

Edge::Edge(int v1, int v2)
{
    this->vertices[0] = v1;
    this->vertices[1] = v2;
}

int Edge::EV(int pos)
{
    return this->vertices[pos];
}
