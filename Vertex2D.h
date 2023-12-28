#ifndef VERTEX2D_H
#define VERTEX2D_H

/**
 * @brief The Vertex2D class represent a point in 2D space
 */

#include <glm/glm.hpp>

class Vertex
{
public:
    Vertex() = default;
    virtual ~Vertex(){}

    void SetVTStar(int pV) {m_VTstar = pV;}
    int VTstar() {return m_VTstar;}
    void VTstar(int pT) {m_VTstar = pT;}

private:
    int m_VTstar;
};


class Vertex2D : public Vertex
{
public:
    ///A constructor method
    Vertex2D();
    ///A constructor method
    Vertex2D(const Vertex2D& pOrig);
    ///A constructor method
    /*!
     * \param x a float argument, representing the x coordinate
     * \param y a float argument, representing the y coordinate
     */
    Vertex2D(double pX, double pY);

    Vertex2D(const glm::vec2& pCoords);
    ///A destructor method
    virtual ~Vertex2D();
    ///
    friend bool operator== (const Vertex2D& p, const Vertex2D &q);
    ///
    friend bool operator!= (const Vertex2D& p, const Vertex2D &q);
    ///A public method that returns the x coordinate
    /*!
     * \return a double value, representing the x coordinate
     */
    double getX();
    double getX() const;

    ///A public method that returns the y coordinate
    /*!
     * \return a double value, representing the y coordinate
     */
    double getY();
    double getY() const;

    ///A public method returning the vertex coordinates
    /*!
     * \return a glm::vec2 with the x and y coordinate of the vertex
     */
    glm::vec2 Coordinates();

    ///A public method that sets the vertex coordinates
    /*!
     * \param pCoords a glm::vec2 containing the coordinates to assign to the vertex
     */
    void SetCoordinates(const glm::vec2& pCoords);
    ///A public method that sets the vertex coordinates
    /*!
     * \param pX a double containing the vertex x coordinate
     * \param pY a double containing the vertex y coordinate
     */
    void SetCoordinates(double pX, double pY);


    ///A public method that sets the x coordinate
    /*!
     * \param pX a double argument, represents the value of the pX coordinate to set
     */
    void setX(double pX);
    ///A public method that sets the x coordinate
    /*!
     * \param pY a double argument, represents the value of the pY coordinate to set
     */
    void setY(double pY);

protected:
    glm::vec2 m_Coordinates;

};

#endif // VERTEX2D_H
