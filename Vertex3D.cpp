#include "Vertex3D.h"

/**
 * @brief Vertex3D::Vertex3D class constructor
 */
Vertex3D::Vertex3D() {
    this->x=0;
    this->y=0;
    this->z=0;
    this->vertexSaliency=0;
    this->vtstar=-1;
}

/**
 * @brief Vertex3D::Vertex3D class construction
 * @param orig
 */
Vertex3D::Vertex3D(const Vertex3D& orig) : Vertex2D(orig) {
    this->z=orig.z;
    //this->vertexSaliency=orig.saliency();
}

Vertex3D::~Vertex3D() {
}

/**
 * @brief Vertex3D::Vertex3D class construction
 * @param x x coordinate
 * @param y y coordinate
 * @param z z coordinate
 */
Vertex3D::Vertex3D(double x, double y, double z){
   this->x = x;
   this->y = y;
   this->z = z;
   this->vtstar = -1;
   this->vertexSaliency=0;
}

void Vertex3D::setZ(double z){
    this->z = z;
}

double Vertex3D::getZ(){
    return this->z;
}

bool operator== (const Vertex3D &p, const Vertex3D &q) {
        return ((p.x == q.x) && (p.y == q.y) && (p.z == q.z));
}

bool operator !=(const Vertex3D& p, const Vertex3D& q) {
        return !(p == q);
}
