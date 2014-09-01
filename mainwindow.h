#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QDir>
#include <QTextEdit>
#include <QLabel>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QFileDialog>

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    QPushButton *ChooseMesh;
    QPushButton *StartSeg;
    QPushButton *Visualize;
    QPushButton *SaveSeg;
    QPushButton *SaveImg;

    QLabel *LabAlpha;
    QLabel *LabNR;
    QLabel *NameFile;
    QLabel *LabIters;
    QLabel *labEta;

    QTextEdit *textAlpha;
    QTextEdit *textNR;
    QTextEdit *textNF;
    QTextEdit *textIters;
    QTextEdit *textEta;
};

#endif // MAINWINDOW_H
