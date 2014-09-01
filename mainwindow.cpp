#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    ChooseMesh = new QPushButton("Choose Mesh", this);
    ChooseMesh->setGeometry(QRect(QPoint(10, 10), QSize(200, 50)));
    StartSeg = new QPushButton("START", this);
    StartSeg->setGeometry(QRect(QPoint(10, 190), QSize(200, 50)));
    Visualize = new QPushButton("Show", this);
    Visualize->setGeometry(QRect(QPoint(380, 10), QSize(200, 50)));
    SaveSeg = new QPushButton("Save segmentation", this);
    SaveSeg->setGeometry(QRect(QPoint(110, 250), QSize(100, 50)));
    SaveImg = new QPushButton("Save Image", this);
    SaveImg->setGeometry(QRect(QPoint(380, 70), QSize(200, 50)));

    LabNR = new QLabel(this);
    LabNR->setGeometry(10, 70, 80, 30);
    LabNR->setText("Regions");
    LabAlpha = new QLabel(this);
    LabAlpha->setGeometry(10, 100, 80, 30);
    LabAlpha->setText("Alpha");
    labEta = new QLabel(this);
    labEta->setGeometry(10, 130, 80, 30);
    labEta->setText("Eta");
    LabIters = new QLabel(this);
    LabIters->setGeometry(10, 160, 80, 30);
    LabIters->setText("Iterations");

    textNR = new QTextEdit(this);
    textNR->setGeometry(100, 70, 90, 30);
    textAlpha = new QTextEdit(this);
    textAlpha->setGeometry(100, 100, 90, 30);
    textNF = new QTextEdit(this);
    textNF->setGeometry(100, 130, 90, 30);
    textIters = new QTextEdit(this);
    textIters->setGeometry(100, 160, 90, 30);


}

MainWindow::~MainWindow()
{
    
}
