#include "Vertex2D.h"

/**
 * @brief Vertex2D::Vertex2D class constructor
 */
Vertex2D::Vertex2D() {
    this->x = 0;
    this->y = 0;
    this->vtstar = -1;
}

/**
 * @brief Vertex2D::Vertex2D class constructor
 * @param orig origin
 */
Vertex2D::Vertex2D(const Vertex2D& orig) {
    this->x = orig.x;
    this->y = orig.y;
    this->vtstar = orig.vtstar;
}

/**
 * @brief Vertex2D::Vertex2D class constructor
 * @param x x coordinate
 * @param y y coordinate
 */
Vertex2D::Vertex2D(double x, double y){
    this->x = x;
    this->y = y;
    this->vtstar = -1;
}

Vertex2D::~Vertex2D() {
}

double Vertex2D::getX(){
    return this->x;
}

double Vertex2D::getY(){
    return this->y;
}


void Vertex2D::setX(double x){
    this->x = x;
}

void Vertex2D::setY(double y){
    this->y = y;
}

int Vertex2D::VTstar()
{
    return this->vtstar;
}

void Vertex2D::VTstar(int vtstar)
{
    this->vtstar = vtstar;
}

bool operator== (const Vertex2D &p, const Vertex2D &q) {
        return ((p.x == q.x) && (p.y == q.y));
}

bool operator !=(const Vertex2D& p, const Vertex2D& q) {
        return !(p == q);
}
