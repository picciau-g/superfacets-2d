/*
 *
 *   2014
 *   Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
 *   Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
 *                 Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
 *
 *   Title:          Fast and scalable mesh superfacets
 *   Submission to Pacific Graphics 2014
 *
 *
 **/

#include "meshvisualizer.h"

MeshVisualizer::MeshVisualizer(QWidget *parent) :
    QWidget(parent)
{
}

/**
 * @brief MeshVisualizer::MeshVisualizer class constructor
 * @param fname name of the mesh file
 * @param sname nae of the segmentation file
 */
MeshVisualizer::MeshVisualizer(string fname, string sname){

    this->meshFile = fname;
    this->segFile = sname;

    initParams();

    findMExtension(fname);
    cout<<"Loaded Mesh "<<fname.c_str()<<endl;
    typeV = 0; //Default normal visualization
    InitStructures();
}

/**
 * @brief MeshVisualizer::InitStructures reserve space for the structures storing information about the object
 */
void MeshVisualizer::InitStructures(){

    /// Initialize the normals
    normals.reserve(mesh.getTopSimplexesNum()*sizeof(Normals));

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        Triangle T = mesh.getTopSimplex(ii);
        Normals n = Normals(mesh.getVertex(T.TV(0)), mesh.getVertex(T.TV(1)),mesh.getVertex(T.TV(2)));
        normals.push_back(n);
    }

    /// Initialize the array with the indices
    clusterIndex = new int[mesh.getTopSimplexesNum()];

    /// Reads the segmentations
    readSegmFromFile(segFile);
    /// Set number of regions
    setNclus();
    /// Mesh centroid
    getCenterMesh();
    /// Bounding box diagonal
    getBBDiagonal();
}

/**
 * @brief MeshVisualizer::findMExtension extract the extension of the mesh file and loads the mesh
 * @param s: mesh name
 */
void MeshVisualizer::findMExtension(string s){

    mesh=Mesh<Vertex3D, Triangle>();
    QString qmeshName=QString::fromStdString(s);
    QStringList ql=qmeshName.split(".");
    if(ql.back().compare("tri")==0)
        Reader::readMeshFile(mesh, s);
    else if(ql.back().compare("off")==0)
        Reader::readOFFMesh(mesh,s);
    mesh.build();
}


/**
 * @brief MeshVisualizer::readSegmFromFile reads the indices of each triangle from file
 * @param segfilename name of the file which stores the segmentation
 */
void MeshVisualizer::readSegmFromFile(string segfilename){

    ifstream ifs;
    string line;

    ifs.open(segfilename.c_str());

    if(ifs.is_open()){

        int indexF = 0;
        while(getline(ifs, line)){

            if(line[0] != '#' && strcmp(line.c_str(), "")){
                unsigned int cluster;
                cluster = atoi(line.c_str());
                clusterIndex[indexF++] = cluster;
            }
        }
    }
}

/**
 * @brief MeshVisualizer::saveImage opens a dialog window to save the segmentation as image
 */
void MeshVisualizer::saveImage(){

    // dialogS *SD = new dialogS;
    // SD->resize(400, 200);

    // QImage QI = this->grab();//grabFrameBuffer();
    // SD->outIm = QI;
    // SD->show();
}

/**
 * @brief MeshVisualizer::getBBDiagonal retrieves lower left and upper right corner of the mesh bounding box
 *          and the length of its diagonal
 */
void MeshVisualizer::getBBDiagonal(){

    // something strange happened if we used infinity
    minCoord.setX(mesh.getVertex(0).getX());
    minCoord.setY(mesh.getVertex(0).getY());
    minCoord.setZ(mesh.getVertex(0).getZ());

    maxCoord.setX(mesh.getVertex(0).getX());
    maxCoord.setY(mesh.getVertex(0).getY());
    maxCoord.setZ(mesh.getVertex(0).getZ());

    for(int ii=1;ii<mesh.getNumVertex();ii++){

        Vertex3D vv = mesh.getVertex(ii);

        //Minimum
        if(vv.getX() < minCoord.getX())
            minCoord.setX(vv.getX());
        if(vv.getY() < minCoord.getY())
            minCoord.setY(vv.getY());
        if(vv.getZ() < minCoord.getZ())
            minCoord.setZ(vv.getZ());

        //Maximum
        if(vv.getX() > maxCoord.getX())
            maxCoord.setX(vv.getX());
        if(vv.getY() > maxCoord.getY())
            maxCoord.setY(vv.getY());
        if(vv.getZ() > maxCoord.getZ())
            maxCoord.setZ(vv.getZ());

        //Retrieve the diagonal
        diagonal = maxCoord.distance(minCoord);
    }
}


/**
 * @brief MeshVisualizer::getCenterMesh retrieves the barycenter of the mesh
 */
void MeshVisualizer::getCenterMesh(){

    float sumX=0.0, sumY=0.0, sumZ=0.0;

    for(int ii=0;ii<mesh.getNumVertex();ii++){

        // sum all vertices coordinates
        sumX += mesh.getVertex(ii).getX();
        sumY += mesh.getVertex(ii).getY();
        sumZ += mesh.getVertex(ii).getZ();

        // divide by number of vertices
        baryCenter.setX(sumX/mesh.getNumVertex());
        baryCenter.setY(sumY/mesh.getNumVertex());
        baryCenter.setZ(sumZ/mesh.getNumVertex());
    }
}

/**
 * @brief MeshVisualizer::setColor decides the current color
 * @param indexT: index of the triangle currently being drawn
 */
void MeshVisualizer::setColor(int indexT){

    ColorMapper CM;
    CM=ColorMapper(numClusters);

    CM.setRed(clusterIndex[indexT]);
    CM.setGreen(clusterIndex[indexT]);
    CM.setBlue(clusterIndex[indexT]);

    glColor3f(CM.red, CM.green, CM.blue);
}

/**
 * @brief MeshVisualizer::draw_triangle draws a triangle of the mesh
 * @param T the triangle to be drawn
 * @param f index of the triangle to be drawn (useful to set the normals)
 */
void MeshVisualizer::draw_triangle(Triangle T, int f){

    int indverts[3];
    for(int i=0;i<3;i++){
        indverts[i]=T.TV(i);
    }
    Vertex3D verts[3];
    for(int i=0;i<3;i++){
        verts[i]=mesh.getVertex(indverts[i]);
    }

    switch(typeV){
    case 0:
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        break;
    case 1:
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        break;
    case 2:
        glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
        glPointSize(3);
        break;
    default:
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }


    glBegin(GL_TRIANGLES);
    for(int ii=0;ii<3;ii++){
        glNormal3f(normals.at(f).getNx(), normals.at(f).getNy(), normals.at(f).getNz());
        glVertex3f(verts[ii].getX(), verts[ii].getY(), verts[ii].getZ());
    }
    glEnd();
}

/**
 * @brief MeshVisualizer::draw_mesh draws the object to screen
 */
void MeshVisualizer::draw_mesh(){

    glMatrixMode(GL_MODELVIEW);
    glColor3f(1.0,1.0,1.0);
    glPushMatrix();

    //positioning the camera
    glRotatef(rotationX, 1.0,0.0,0.0);
    glRotatef(rotationY, 0.0,1.0,0.0);
    glRotatef(rotationZ, 0.0,0.0,1.0);

    glPopMatrix();

    //position the object
    glTranslatef(-objPX, -objPY, -objPZ);
    glRotatef(objRZ, 0.0, 0.0,1.0);
    glRotatef(objRY, 0.0,1.0,0.0);
    glRotatef(objRX, 1.0,0.0,0.0);

    for(faceind ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        //Display the ii-th triangle
        Triangle T = mesh.getTopSimplex(ii);
        setColor(ii);
        draw_triangle(T, ii);

        /// If visualization mode is not points, mark edges between clusters
        if(typeV!=2){
            for(int jj=0;jj<3;jj++){

                int neigh = T.TT(jj);

                if(clusterIndex[ii] != clusterIndex[neigh]){

                    Edge* ee = T.TE(jj);
                    Vertex3D v0 = mesh.getVertex(ee->EV(0));
                    Vertex3D v1 = mesh.getVertex(ee->EV(1));

                    glLineWidth(2.0);
                    glColor3f(0.0,0.0,0.0);

                    glBegin(GL_LINES);
                    glVertex3f(v0.getX(), v0.getY(), v0.getZ());
                    glVertex3f(v1.getX(), v1.getY(), v1.getZ());
                    glEnd();
                }
            }
        }
    }
}

/**
 * @brief MeshVisualizer::initializeGL opengl initialization function
 */
void MeshVisualizer::initializeGL(){

    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_TEXTURE_2D);
    glClearColor(1.0,1.0,1.0,1.0);

}

/**
 * @brief MeshVisualizer::resizeGL opengl window function
 * @param width the width of the window
 * @param height the height of the window
 */
void MeshVisualizer::resizeGL(int width, int height){
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    GLfloat x = GLfloat(width) / GLfloat(height);
    gluPerspective(25.0, x, 0.2, 200);
    glMatrixMode(GL_MODELVIEW);
}

/**
 * @brief MeshVisualizer::paintGL draws the scene
 */
void MeshVisualizer::paintGL(){

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);

    glLoadIdentity();

    /// Observer is here
    Vertex3D eye;
    eye.setX(baryCenter.getX() + (maxCoord.getX() - baryCenter.getX())*3.5*sqrt(3));
    eye.setY(baryCenter.getY() + (maxCoord.getY() - baryCenter.getY())*3.5);
    eye.setZ(baryCenter.getZ() + (maxCoord.getZ() - baryCenter.getZ())*3.5);

    gluLookAt(eye.getX(), eye.getY(), eye.getZ(), baryCenter.getX(), baryCenter.getY(), baryCenter.getZ(), 0, 1, 0);

    /// draw the object
    draw_mesh();
}

/**
 * @brief MeshVisualizer::mousePressEvent mouse click manager
 * @param e event
 */
void MeshVisualizer::mousePressEvent(QMouseEvent *e){

    lPos = e->pos();
}

/**
 * @brief MeshVisualizer::mouseMoveEvent manage mouse move event
 * @param e
 */
void MeshVisualizer::mouseMoveEvent(QMouseEvent *e){

    //Amplitude of the move
    GLfloat dx = GLfloat((e->x() - lPos.x()))/width();
    GLfloat dy = GLfloat((e->y() - lPos.y()))/height();

    if(e->buttons() & Qt::LeftButton){
        objRX += 180*dx;
        objRY += 180*dy;
    }
    else if(e->buttons() & Qt::RightButton){
        objRX += 180*dy;
        objRZ += 180*dx;
    }
    update();
    lPos = e->pos();
}

/**
 * @brief MeshVisualizer::keyPressEvent keyboard interaction
 * @param e
 */
void MeshVisualizer::keyPressEvent(QKeyEvent *e){

    switch(e->key()){

    //close window
    case Qt::Key_Escape:
    case Qt::Key_Q:
        close();
        break;

    case Qt::Key_Up:
        objPY+=0.05*diagonal;
        break;

    case Qt::Key_Down:
        objPY-=0.05*diagonal;
        break;

    case Qt::Key_Left:
        objPZ+=0.05*diagonal;
        break;

    case Qt::Key_Right:
        objPZ-=0.05*diagonal;
        break;

    case Qt::Key_E:
        objPX -= 0.05*diagonal;
        break;

    case Qt::Key_X:
        objPX += 0.05*diagonal;
        break;

    case Qt::Key_S:
        saveImage();
        break;

    case Qt::Key_F:         /// Triangle visualization (default)
        typeV=0;
        break;
    case Qt::Key_L:         /// Line visualization
        typeV=1;
        break;
    case Qt::Key_P:         /// Point visualization
        typeV=2;
        break;

    default:
        e->ignore();
        break;
    }
    update();
}
