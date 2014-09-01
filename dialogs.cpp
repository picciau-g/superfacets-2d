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

#include "dialogs.h"

/**
 * @brief dialogS::dialogS creates the dialog window
 * @param parent
 */
dialogS::dialogS(QWidget *parent) :
    QWidget(parent)
{
    nameToGive = new QLabel(this);
    nameToGive->setText("Save as:");
    nameToGive->setGeometry(10, 20, 100, 30);

    imgName = new QTextEdit(this);
    imgName->setGeometry(150, 20, 200, 30);

    bSave = new QPushButton("Save", this);
    bSave->setGeometry(150, 100, 100, 50);

    connect(imgName, SIGNAL(textChanged()), this, SLOT(nameGiven()));
    connect(bSave, SIGNAL(clicked()), this, SLOT(buttonSaveClicked()));
}

/**
 * @brief dialogS::nameGiven set the name of the image to be saved
 */
void dialogS::nameGiven(){

    this->nameI = imgName->toPlainText();
}

/**
 * @brief dialogS::buttonSaveClicked actually saves the content of the visualizer to a file
 */
void dialogS::buttonSaveClicked(){

    outIm.save(imgName->toPlainText());
    this->close();
}

