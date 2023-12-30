#include "Vertex3D.h"

/**
 * @brief Vertex3D::Vertex3D class constructor
 */
Vertex3D::Vertex3D()
    :
    m_Coordinates(0.0,0.0,0.0)
    , vertexSaliency(0.0)
    , m_VTstar(-1)
{

}

/**
 * @brief Vertex3D::Vertex3D class construction
 * @param pOrig
 */
Vertex3D::Vertex3D(const Vertex3D& pOrig)
    :
    m_Coordinates(pOrig.m_Coordinates)
    , vertexSaliency(pOrig.vertexSaliency)
    , m_VTstar(pOrig.m_VTstar)
{

}

Vertex3D::~Vertex3D() {
}

/**
 * @brief Vertex3D::Vertex3D class construction
 * @param pX pX coordinate
 * @param pY pY coordinate
 * @param pZ pZ coordinate
 */
Vertex3D::Vertex3D(double pX, double pY, double pZ)
:
    m_Coordinates(pX, pY, pZ)
    , vertexSaliency(0.0)
    , m_VTstar(-1)
{

}


Vertex3D::Vertex3D(const glm::vec3& pCoords)
    :
    m_Coordinates(pCoords)
    , vertexSaliency(0.0)
    , m_VTstar(-1)
{

}


void Vertex3D::setX(double pX)
{
    this->m_Coordinates.x = pX;
}

void Vertex3D::setY(double pY)
{
    this->m_Coordinates.y = pY;
}

void Vertex3D::setZ(double pZ)
{
    this->m_Coordinates.z = pZ;
}


double Vertex3D::getX()
{
    return m_Coordinates.x;
}

double Vertex3D::getX() const
{
    return m_Coordinates.x;
}


double Vertex3D::getY()
{
    return m_Coordinates.y;
}

double Vertex3D::getY() const
{
    return m_Coordinates.y;
}

double Vertex3D::getZ()
{
    return m_Coordinates.z;
}

double Vertex3D::getZ() const
{
    return m_Coordinates.z;
}


glm::vec3 Vertex3D::GetCoordinates() const
{
    return m_Coordinates;
}


double Vertex3D::SquaredMagnitude()
{
    return glm::length2(m_Coordinates);
}

double Vertex3D::SquaredDistance(const Vertex3D& pOth)
{
    return glm::length2(pOth.GetCoordinates() - m_Coordinates);
}


double Vertex3D::prodscal(const Vertex3D& v1,const Vertex3D& v2)
{
    return glm::dot(v1.GetCoordinates(), v2.GetCoordinates());
}


double Vertex3D::prodscal(const Vertex3D& pOther)
{
    return glm::dot(m_Coordinates, pOther.GetCoordinates());
}


bool operator== (const Vertex3D &p, const Vertex3D &q) {
    return p.GetCoordinates() == q.GetCoordinates();
}

bool operator !=(const Vertex3D& p, const Vertex3D& q) {
        return !(p == q);
}
