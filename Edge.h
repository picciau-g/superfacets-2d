#ifndef EDGE_H
#define EDGE_H

/**
 * @brief The Edge class stores minimal information about edges of a mesh
 */

class Edge
{
public:
    Edge();
    Edge(int v1, int v2);

    int EV(int pos);

    inline bool operator== (const Edge &p)
    {
        bool b[2];
        b[0] = false; b[1] = false;
        for(int i=0;i<2;i++)
        {
            for(int j=0;j<2;j++)
            {
                if(!b[j] && p.vertices[i]==this->vertices[j])
                {
                    b[j] = true;
                    break;
                }
            }
        }

        return b[0] && b[1];
    }

    inline bool operator!= (Edge &p)
    {
       return !(p==(*this));
    }

private:
    int vertices[2];
};


#endif // EDGE_H
