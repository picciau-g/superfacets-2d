#ifndef VERTEX2D_H
#define VERTEX2D_H

/**
 * @brief The Vertex2D class represent a point in 2D space
 */

class Vertex2D
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
    ///
    int VTstar();
    ///
    void VTstar(int vtstar);
protected:
    ///A protected variable representing the m_X coordinate of the point
    double m_X;
    ///A protected variable representing the m_Y coordinate of the point
    double m_Y;
    ///
    int m_VTstar;
};

#endif // VERTEX2D_H
