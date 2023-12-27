#include "Vertex3D.h"

/**
 * @brief Vertex3D::Vertex3D class constructor
 */
Vertex3D::Vertex3D() :
    Vertex2D(0,0)
    , m_Z(0)
    , vertexSaliency(0.0)
{

}

/**
 * @brief Vertex3D::Vertex3D class construction
 * @param pOrig
 */
Vertex3D::Vertex3D(const Vertex3D& pOrig) : Vertex2D(pOrig)
{
    this->m_Z=pOrig.m_Z;
    //this->vertexSaliency=orig.saliency();
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
    Vertex2D(pX,pY)
    , m_Z(pZ)
    , vertexSaliency(0.0)
{

}

void Vertex3D::setZ(double pZ)
{
    this->m_Z = pZ;
}

double Vertex3D::getZ()
{
    return m_Z;
}

double Vertex3D::getZ() const
{
    return m_Z;
}


double Vertex3D::Norm()
{
    return sqrt(m_X*m_X + m_Y*m_Y + m_Z*m_Z);
}

double Vertex3D::Norm(const Vertex3D& pOth)
{
    return sqrt((pOth.getX() - m_X)*(pOth.getX() - m_X) +
                + (pOth.getY() - m_Y)*(pOth.getY() - m_Y) +
                + (pOth.getZ() - m_Z)*(pOth.getZ() - m_Z));
}

double Vertex3D::prodscal(const Vertex3D& pOther)
{
    return ((pOther.getX())*(m_X))+((pOther.getY())*(m_Y))+((pOther.getZ())*(m_Z));
}


bool operator== (const Vertex3D &p, const Vertex3D &q) {
    return ((p.m_X == q.m_X) && (p.m_Y == q.m_Y) && (p.m_Z == q.m_Z));
}

bool operator !=(const Vertex3D& p, const Vertex3D& q) {
        return !(p == q);
}
